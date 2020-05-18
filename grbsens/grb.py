import numpy as np


class grb:
    """This is a class to store all the information about a GRB that is to be run

    :param num_jobs: Total number of jobs, defaults to 1
    :type num_jobs: (int, NoneType), optional
    :param total_time: Total duration (s) over which to run sims, defaults to None
    :type total_time: (float, NoneType), optional
    :param delta_t: Duration of time steps (s) over which to iterate sims, defaults to 1
    :type delta_t: (float, NoneType), optional
    :param lower_limit: lower energy limit in TeV, defaults to 0.02 TeV
    :type lower_limit: float, optional
    :param upper_limit: upper energy limit in Tev, defaults to 10 TeV
    :type upper_limit: float, optional
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
    """

    def __init__(
            self,
            num_jobs=1,
            total_time=None,
            delta_t=1.,
            lower_limit=0.02,
            upper_limit=10,
            bins=1,
            irf="North_0.5h",
            init_time=1.,
            sigma=5.,
            offset=0.,
    ):
        """Constructor method
        """

        # initialize parameters dictionary
        self.params = {
            "num_jobs": num_jobs,
            "total_time": total_time,
            "delta_t": delta_t,
            "lower_limit": lower_limit,
            "upper_limit": upper_limit,
            "bins": bins,
            "irf": irf,
            "init_time": init_time,
            "sigma": sigma,
            "offset": offset,
        }

        # check inputs
        self._check_inputs()

        # get time steps
        self._get_time_steps()

    def _check_inputs(self):
        """Check the validity of the inputs upon class initialization"""

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

        # add stop time to params
        self.params["stop_time"] = stop_time

        print(f"Running from t0={start_time}s to t1={stop_time}s "
              f"for a total duration of t={self.params['total_time']} "
              f"with {self.params['num_jobs']} time steps of dt={time_step}s each")

        self.times = times


if __name__ == "__main__":
    my_grb = grb(init_time=0, total_time=25, delta_t=0.5, num_jobs=None)

    print("done")
