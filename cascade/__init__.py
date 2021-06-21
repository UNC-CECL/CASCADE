import pkg_resources

from .cascade import Cascade

__version__ = pkg_resources.get_distribution("cascade").version
__all__ = ["Cascade"]

del pkg_resources
