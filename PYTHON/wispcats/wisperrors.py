

class Error(Exception):
    """Base class for exceptions raised in wispcats."""
    pass


class WISPError(Error):
    """ """
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


