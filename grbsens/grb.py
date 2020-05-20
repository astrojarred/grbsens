import os
import numpy as np
import pandas as pd

import gammalib
import ctools
import cscripts


class grb:
    """This is a class to store all the information about a GRB that is to be run

    :param num_jobs: Total number of jobs, defaults to 1
    :type num_jobs: (int, NoneType), optional
    :param total_time: Total duration (s) over which to run sims, defaults to None
    :type total_time: (float, NoneType), optional
    :param delta_t: Duration of time steps (s) over which to iterate sims, defaults to 1
    :type delta_t: (float, NoneType), optional
    :param emin: lower energy limit in TeV, defaults to 0.02 TeV
    :type emin: float, optional
    :param emax: upper energy limit in Tev, defaults to 10 TeV
    :type emax: float, optional
    :param bins: number of energy bins, defaults to 1
    :type bins: int, optional
    :param irf: input ctools IRF, defaults to "North_0.5h"
    :type irf, string, optional
    :param init_time: initial observation time, defaults to 1
    :type init_time: float, optional
    :param sigma: significance, defaults to 5 sigma
    :type sigma: float, optional
    :param offset: offset from the pointing coordinates
    :type offset: float, optional
    TODO: add rest of the parameters
    """

    def __init__(
            self,
            input_model,
            num_jobs=1,
            total_time=None,
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
        """Constructor method
        """

        self.input_model = input_model

        # initialize parameters dictionary
        self.params = {
            "num_jobs": num_jobs,
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

        """initialize results dict:
        results = {
            index = {
                job number  : xxx
                output file : yyy.txt
                log file    : zzz.log
            }
        """
        self.results = {}

        # check inputs
        self._check_inputs()

        # get time steps
        self._get_time_steps()

    def _check_inputs(self):
        """Check the validity of the inputs upon class initialization"""

        # check sensitivity type is either integral or differential
        sens_type = self.params["sens_type"].lower()
        if sens_type != "integral" and sens_type != "differential":
            raise AttributeError("Parameter `sens_type` must be"
                                 " either 'Integral' or 'Differential'.")
        # capitalize sens_type
        self.params["sens_type"] = sens_type.capitalize()

        # check delta_t, num_jobs, total time
        number_of_params = 0
        param_not_set = ""
        for param in ["num_jobs", "total_time", "delta_t"]:
            if self.params[param] is not None:
                number_of_params += 1
            else:
                param_not_set = param

        # confirm there are only 2/3 parameters given
        if number_of_params > 2:
            raise AttributeError("Please do not provide more than 2 of the 3 parameters: "
                                 "`num_jobs`, `total_time`, `delta_t`")

        # calculate third parameter
        if param_not_set == "num_jobs":
            self.params["num_jobs"] = int(self.params["total_time"] / self.params["delta_t"])
        elif param_not_set == "total_time":
            self.params["total_time"] = self.params["num_jobs"] * self.params["delta_t"]
        elif param_not_set == "delta_t":
            self.params["delta_t"] = self.params["total_time"] / self.params["num_jobs"]

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
        outfile = f"./outputs/sensi-{self.params['sigma']}sigma_obstime-{duration}_irf-{self.params['irf']}.txt"
        logfile = f"./logs/sensi-{self.params['sigma']}sigma_obstime-{duration}_irf-{self.params['irf']}.log"

        # load input model
        models = gammalib.GModels(self.input_model)
        models.save(self.input_model)
        sen["inmodel"] = self.input_model

        sen["duration"] = duration
        sen["outfile"] = outfile
        sen["logfile"] = logfile

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

        # add duration as a column
        results['duration'] = [duration]

        # add outputs to dicts
        self.results[job_number] = dict(
            duration=duration,
            table=results,
            log=logfile,
        )

    def execute(self):

        for job_number, duration in enumerate(self.times):
            self._calculate_sensitivity(job_number=job_number, duration=duration)
            print(f"Done with duration={duration}s\n")


if __name__ == "__main__":
    my_grb = grb(init_time=0, total_time=4, delta_t=1, num_jobs=None)

    print("done")
