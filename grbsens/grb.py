class grb:
    """This is a class to store all the information about a GRB that is to be run

    :param num_runs: Total number of runs, defaults to 1
    :type num_runs: int, optional
    :param num_jobs: Total number of jobs, defaults to 1
    :type num_jobs: int, optional
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
        num_runs=1,
        num_jobs=1,
        lower_limit=0.02,
        upper_limit=10,
        bins=1,
        irf="North_0.5h",
        init_time=1,
        delta_t=1,
        sigma=5,
        offset=0.0,
    ):
        """Constructor method
        """
        # initialize parameters dictionary
        self.params = {
            "num_runs": num_runs,
            "num_jobs": num_jobs,
            "lower_limit": lower_limit,
            "upper_limit": upper_limit,
            "bins": bins,
            "irf": irf,
            "init_time": init_time,
            "delta_t": delta_t,
            "sigma": sigma,
            "offset": offset,
        }

