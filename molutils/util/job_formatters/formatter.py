class Formatter(object):
    def format(self, type, method, guess_charge=False):
        """
        This method is a wrapper function to call the appropriate formatter method
        based on the type of calculation and method of calculation.
        :param type: calculation type (e.g. energy)
        :param method: calculation method (e.g. MP2)
        :return: the formatted result
        """
        raise NotImplemented("All formatters must implement the format method")