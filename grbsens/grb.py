import os
import numpy as np
import pandas as pd

import gammalib
import cscripts


class grb:
    """This is a class to store all the information about a GRB that is to be run."""

    def __init__(
            self,
            input_model,
            total_time=1.,
            delta_t=1.,
            emin=0.03,
            emax=10,
            bins=1,
            irf="North_0.5h",
            init_time=1.,
            sigma=5.,
            offset=0.,
            binsz=0.2,
            sens_type="Integral",
            rad=2.25,
            caldb="prod2",
            src_name="GRB",
    ):

        self.input_model = input_model

        # initialize parameters dictionary
        self.params = {
            "total_time": total_time,
            "delta_t": delta_t,
            "emin": emin,
            "emax": emax,
            "bins": bins,
            "irf": irf,
            "init_time": init_time,
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

        # check inputs
        self._check_inputs()

        # get time steps
        self._get_time_steps()

    def __len__(self):
        """Returns number of time steps"""
        return len(self.times)

    def __getitem__(self, item):
        """Return parameter"""
        return self.params[item]

    def __setitem__(self, key, value):
        """Set a parameter manually"""
        self.params[key] = value

    def _check_inputs(self):
        """Check the validity of the inputs upon class initialization"""

        # check sensitivity type is either integral or differential
        sens_type = self.params["sens_type"].lower()
        if sens_type != "integral" and sens_type != "differential":
            raise AttributeError("Parameter `sens_type` must be"
                                 " either 'Integral' or 'Differential'.")
        # capitalize sens_type
        self.params["sens_type"] = sens_type.capitalize()

        # calculate the number of jobs parameter
        self.params["num_jobs"] = int(self.params["total_time"] / self.params["delta_t"])

    def _get_time_steps(self):
        """Create the time steps for the class."""

        start_time = self.params["init_time"]
        stop_time = self.params["init_time"] + self.params["total_time"]
        time_step = self.params["delta_t"]
        times = np.arange(start_time + time_step, stop_time + time_step, time_step)

        # add stop time to params2
        self.params["stop_time"] = stop_time

        print(f"Running from t0={start_time}s to t1={stop_time}s "
              f"for a total duration of t={self.params['total_time']} "
              f"with {self.params['num_jobs']} time steps of dt={time_step}s each")

        self.times = times

    def _calculate_sensitivity(self, job_number, duration):
        """Run the `cssens` ctools module based on the given input"""

        # set duration to a float
        duration = float(duration)

        print(f"Running `cssens` job #{job_number} for "
              f"{self.params['src_name']} for a duration of {duration}s")

        # create cssens object
        sen = cscripts.cssens()

        # calculate outfile and logfile names
        cwd = os.path.abspath('')  # current working directory
        outfile = f"{cwd}/outputs/sensi-{self.params['sigma']}sigma_obstime-{duration}_irf-{self.params['irf']}.txt"
        logfile = f"{cwd}/logs/sensi-{self.params['sigma']}sigma_obstime-{duration}_irf-{self.params['irf']}.log"

        # load input model
        models = gammalib.GModels(self.input_model)
        models.save(self.input_model)
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

        # run cssens
        sen.execute()

        # import results into a dataframe
        results = pd.read_csv(outfile)

        # add duration and job number as a column
        results['duration'] = [duration]
        results['job_number'] = [job_number]

        # add output and log filepaths as columns
        results['output_file'] = [outfile]
        results['log_file'] = [logfile]

        # add output pandas df to results dictionary
        self.results[job_number] = results

    def save_to_csv(self, filepath=None):
        """Save results to a csv"""

        if filepath is None:
            cwd = os.path.abspath('')  # current working directory
            start, stop = min(self.times), max(self.times)  # get start and stop times
            filepath = f"{cwd}/outputs/sensi-{self.params['sigma']}sigma_t{start}s-t{stop}s_irf-{self.params['irf']}.csv"

        # save as csv
        self.output.to_csv(filepath)
        print(f"\nOutput written to {filepath}\n")

    def execute(self, write_to_file=True, output_filepath=None):
        """Run `cssens` once for each job"""

        for job_number, duration in enumerate(self.times):
            self._calculate_sensitivity(job_number=job_number, duration=duration)
            print(f"Done with duration={duration}s\n")

        # concatenate results
        self.output = pd.concat(self.results, ignore_index=True).set_index("job_number")

        # write to csv file
        if write_to_file:
            self.save_to_csv(filepath=output_filepath)


if __name__ == "__main__":
    my_grb = grb(input_model="grb.xml", init_time=0, total_time=4, delta_t=1, num_jobs=None)

    print("done")
