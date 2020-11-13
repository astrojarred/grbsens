import os

import numpy as np
import pandas as pd
import numbers
from matplotlib import pyplot as plt
import multiprocessing as mp
from tqdm.auto import tqdm
from pathlib import Path

import cscripts


class grb:
    """This is a class to store all the information about a GRB that is to be run."""

    def __init__(
            self,
            input_model,
            start_time=0.,
            stop_time=1.,
            delta_t=1.,
            log_steps=None,
            emin=0.03,
            emax=10,
            bins=1,
            irf="North_0.5h",
            npix=100,
            sigma=5.,
            offset=0.,
            binsz=0.2,
            sens_type="Integral",
            rad=2.25,
            caldb="prod2",
            src_name="GRB",
    ):
        """Initialize class. TODO: add parameters."""
        self.input_model = input_model

        # initialize parameters dictionary
        self.params = {
            "stop_time": stop_time,
            "delta_t": delta_t,
            "log_steps": log_steps,
            "emin": emin,
            "emax": emax,
            "bins": bins,
            "irf": irf,
            "npix": npix,
            "start_time": start_time,
            "sigma": sigma,
            "offset": offset,
            "binsz": binsz,
            "sens_type": sens_type,
            "rad": rad,
            "caldb": caldb,
            "src_name": src_name,

        }

        # initialize results dictionary and output
        self.results = {}
        self.output = None
        self.std_output_df = None

        # check inputs
        self._check_inputs()

        # get time steps
        self.times = np.array([])
        self.timeframes = {}
        self._get_time_steps()

    def __len__(self):
        """Returns number of time steps."""
        return len(self.times)

    def __getitem__(self, item):
        """Return parameter."""
        return self.params[item]

    def __setitem__(self, key, value):
        """Set a parameter manually."""
        self.params[key] = value

    def _check_inputs(self):
        """Check the validity of the inputs upon class initialization."""
        # check sensitivity type is either integral or differential
        sens_type = self.params["sens_type"].lower()
        if sens_type != "integral" and sens_type != "differential":
            raise AttributeError("Parameter `sens_type` must be"
                                 " either 'Integral' or 'Differential'.")
        # capitalize sens_type
        self.params["sens_type"] = sens_type.capitalize()

    def _get_time_steps(self):
        """Create the time steps for the class."""

        # catch normal time steps
        if isinstance(self.params["delta_t"], numbers.Number):

            start_time = self.params["start_time"]
            stop_time = self.params["stop_time"]
            time_step = self.params["delta_t"]
            times = np.arange(start_time + time_step, stop_time + time_step, time_step)

            print(f"Running from t0={start_time}s to t1={stop_time}s "
                  f"with time steps of dt={time_step}s each")

            self.times = times

        # catch arrays for manual time steps
        elif isinstance(self.params["delta_t"], (list, np.ndarray)):

            self.times = self.params["delta_t"]
            self.params["start_time"] = min(self.times)
            self.params["stop_time"] = max(self.times)

        elif self.params["delta_t"] == "custom":
            # create time frames dict
            self.timeframes = {}

            print("Add time frames with custom time steps using "
                  "`grb.add_timeframe(start, stop, time_step)` in seconds.")

        elif self.params["delta_t"] == "log":

            if self.params["log_steps"] is None:
                raise AttributeError("For log mode, please specify the initial time, total time, "
                                     "and number of time steps with `start_time`, `stop_time`,"
                                     "and `log_steps`, respectively.")

            start_time = self.params["start_time"]
            stop_time = self.params["stop_time"]
            log_steps = self.params["log_steps"]

            # catch invalid start times
            if np.log10(start_time) == -np.inf:
                raise AttributeError("In log mode, `start_time` cannot be equal to 0.")

            times = np.logspace(np.log10(start_time), np.log10(stop_time), log_steps)

            print(f"Running from t0={start_time}s to t1={stop_time}s "
                  f"with {log_steps} time steps on a log scale.")

            self.times = times

        else:
            raise AttributeError("Choose a either a single value for `delta_t` or select 'custom' to"
                                 "add custom time frames with unique time steps.")

    def add_timeframe(self, start, stop, time_step):

        key = self._get_next_timeframe_key()

        self.timeframes[key] = {
            "start": start,
            "stop": stop,
            "dt": time_step,
        }

        self._get_all_time_steps()
        print(f"Added time frame #{key} from {start}s to {stop}s with time step {time_step}s.")

    def reset_timeframes(self):

        self.times = np.array([])
        self.timeframes = {}

        print("Successfully reset time frames.\n"
              "Add time frames with custom time steps using "
              "`grb.add_timeframe(start, stop, time_step)` in seconds.")

    def _get_next_timeframe_key(self):

        keys = np.fromiter(self.timeframes.keys(), dtype=int)
        if len(keys) == 0:
            return 0
        else:
            return np.max(keys) + 1

    def _get_all_time_steps(self):

        all_times = np.array([])

        for val in self.timeframes.values():
            all_times = np.append(
                all_times,
                np.arange(val["start"] + val["dt"], val["stop"] + val["dt"], val["dt"]),
            )

        self.times = all_times

    def _calculate_sensitivity(self, job_number, duration, cwd=None, parallel_results=None,
                               nthreads=1, load_results=False, verbose=True):
        """Run the `cssens` ctools module based on the given input."""
        # set duration to a float
        duration = float(duration)
        if verbose:
            print(f"Running `cssens` job #{job_number} for "
                  f"{self.params['src_name']} for a duration of {duration}s")

        # create cssens object
        sen = cscripts.cssens()

        outfile = f"{cwd}/cssens_outputs/grbsens-{self.params['sigma']}" \
                  f"sigma_obstime-{duration}_irf-{self.params['irf']}.txt"
        logfile = f"{cwd}/cssens_logs/grbsens-{self.params['sigma']}" \
                  f"sigma_obstime-{duration}_irf-{self.params['irf']}.log"

        # run cssens
        if not load_results or not Path(outfile).is_file():
            # load input model
            sen["inmodel"] = self.input_model

            # set parameters that change each loop
            sen["duration"] = duration
            sen["outfile"] = outfile
            sen["logfile"] = logfile

            # set global parameters
            sen["srcname"] = self.params["src_name"]
            sen["caldb"] = self.params["caldb"]
            sen["irf"] = self.params["irf"]
            sen["rad"] = self.params["rad"]
            sen["emin"] = self.params["emin"]
            sen["emax"] = self.params["emax"]
            sen["type"] = self.params["sens_type"]
            sen["sigma"] = self.params["sigma"]
            sen["bins"] = self.params["bins"]
            sen["binsz"] = self.params["binsz"]
            sen["offset"] = self.params["offset"]
            sen["npix"] = self.params["npix"]

            # set number of cores used for energy bins
            sen["nthreads"] = nthreads
            # set chatter to max
            sen["chatter"] = 4

            sen.execute()

        # format results into pandas DataFrame
        results = self._results_to_df(outfile, duration, job_number, logfile)

        # add output pandas df to results dictionary
        self._save_results(results=results, job_number=job_number, parallel_results=parallel_results)

        if verbose:
            # print success message
            print(f"Done with job #{job_number}, duration={duration}s\n")

    @staticmethod
    def _results_to_df(outfile, duration, job_number, logfile):
        """Format results of _calculate_sensitivity into a pandas DataFrame."""
        # import results into a dataframe
        results = pd.read_csv(outfile)

        # add duration and job number as a column
        results['duration'] = [duration]
        results['job_number'] = [job_number]

        # add output and log file paths as columns
        results['output_file'] = [outfile]
        results['log_file'] = [logfile]

        return results

    def _save_results(self, results, job_number, parallel_results=None):
        """Write results of _calculate sensitivity to the grb class."""
        if parallel_results is None:
            self.results[job_number] = results
        else:
            parallel_results[job_number] = results

    def save_to_csv(self, filename=None, cwd=None):
        """Save results to a csv."""
        if filename is None:
            if cwd is None:
                cwd = os.path.abspath('')  # current working directory
            start, stop = min(self.times), max(self.times)  # get start and stop times
            filename = f"grbsens-{self.params['sigma']}sigma_t{start}s-t{stop}s_irf-{self.params['irf']}.txt"

        filepath = f"{cwd}/grbsens_results/{filename}"

        # The following steps to match the standard output formatting
        # copy output df
        f = self.output.copy(deep=True)

        # Create a new column called Obs time
        f["Obs time"] = f.index

        # Select which columns to print
        f = f[["Obs time", "crab_flux", "photon_flux", "energy_flux", "sensitivity"]]

        # round outputs to 6 decimal places in scientific notation
        for c in ["crab_flux", "photon_flux", "energy_flux", "sensitivity"]:
            f[c] = np.array([f"{i:.6e}" for i in f[c]])

        # space out observations times
        f["Obs time"] = np.array([f"{i:<9}" for i in f["Obs time"]])

        # write comment lines with column names and units
        with open(filepath, 'w') as file:
            file.write("#\n")
            file.write("#Obs time\tcrab_flux\tphoton_flux\tenergy_flux\tsensitivity\n")
            file.write("#s\tcrab units\tph/cm2/s\terg/cm2/s\terg/cm2/s\n")

        # write to csv
        f.to_csv(filepath, sep="\t", index=False, header=False, mode='a')

        # save standard output dataframe source to class
        self.std_output_df = f.copy(deep=True)

        print(f"\nOutput written to {filepath}\n")

    def execute(self, write_to_file=True, output_filename=None, cwd=None, parallel=False,
                ncores=1, nthreads=1, load_results=False, verbose=False):
        """Run `cssens` once for each job."""

        # create output folders if they don't already exist
        if cwd is None:
            cwd = os.path.abspath('')  # set current working directory to execution directory

        Path(f"{cwd}/cssens_outputs").mkdir(parents=True, exist_ok=True)
        Path(f"{cwd}/cssens_logs").mkdir(parents=True, exist_ok=True)
        Path(f"{cwd}/grbsens_results").mkdir(parents=True, exist_ok=True)

        if not parallel:
            for job_number, duration in tqdm(enumerate(self.times),
                                             total=len(self.times),
                                             desc=f'{self.params["src_name"]}'):
                self._calculate_sensitivity(job_number=job_number, duration=duration, cwd=cwd,
                                            nthreads=nthreads, load_results=load_results, verbose=verbose)
        elif parallel:
            # run in parallel with asynchronous pooling
            # check that selected cores is not too many
            if ncores > mp.cpu_count():
                raise AttributeError(f"Selected quantity of cores {ncores} "
                                     f"is greater than available cores {mp.cpu_count()}.")
            # set up pool with ncores CPUs
            pool = mp.Pool(ncores)

            print(f"Running {len(self.times)} jobs in parallel across {ncores} cores:")

            # set up results storage manager
            manager = mp.Manager()
            parallel_results = manager.dict()

            # initialize progress bar
            progress_bar = tqdm(total=len(self.times), desc=f'{self.params["src_name"]}')

            # run loop
            for job_number, duration in enumerate(self.times):
                pool.apply_async(self._calculate_sensitivity,
                                 args=(job_number, duration, cwd, parallel_results, nthreads, load_results, verbose),
                                 callback=lambda _: progress_bar.update(1))

            # Close Pool and let all the processes complete
            pool.close()
            pool.join()  # postpones the execution of next line of code until all processes in the queue are done.

            # save parallel results to class
            self.results = dict(parallel_results)

            # print success message
            print(f"Done running {len(self.times)} jobs in parallel across {ncores} cores!")

        # concatenate results
        self.output = pd.concat(self.results, ignore_index=True).set_index("job_number")

        # set durations to integers
        self.output.duration = [int(i) for i in self.output.duration]

        # set duration as index
        self.output.set_index("duration", inplace=True)

        # sort results by index
        self.output.sort_index(inplace=True)

        # write to csv file
        if write_to_file:
            self.save_to_csv(filename=output_filename, cwd=cwd)

    def plot_results(self, logx=True, logy=True):
        """Plot results on a duration vs sensitivity scatter."""
        plt.figure(figsize=(9, 6))

        # log x and y acxis
        if logy:
            y = np.log10(np.array(self.output.sensitivity))
            y_label = "$\log_{10}$ Sensitivity [erg/cm2/s]"
        else:
            y = np.array(self.output.sensitivity)
            y_label = "Sensitivity [erg/cm2/s]"

        if logx:
            x = np.log10(np.array(self.output.index))
            x_label = "$\log_{10}$ Duration [s]"
        else:
            x = np.array(self.output.index)
            x_label = "Duration [s]"

        plt.scatter(x, y)
        plt.xlabel(x_label, size=15)
        plt.ylabel(y_label, size=15)
        # return fig


if __name__ == "__main__":
    # initialize class
    my_grb = grb(input_model="grb.xml", start_time=0, stop_time=4, delta_t=1)

    # execute grbsens, skip actual running
    my_grb.execute(write_to_file=False, parallel=True, ncores=10)

    print("done")
