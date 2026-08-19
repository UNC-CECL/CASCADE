"""Generic CASCADE grid geometry: real domains padded by buffers on each end.

A CASCADE/BRIE run pads its GIS-numbered real domains with buffer domains
on each side so alongshore diffusion doesn't see an artificial edge at the
real-domain boundary. DomainGeometry's field defaults happen to match the
Hatteras Island hindcast (90 real domains, 15-domain buffers, 500 m
spacing, starting at GIS ID 1) -- override them for a different site or
grid resolution; nothing else in this module assumes Hatteras' numbers.
"""

import dataclasses


@dataclasses.dataclass(frozen=True)
class DomainGeometry:
    """Describes the padded CASCADE domain array.

    Attributes:
        num_real_domains: Count of real, GIS-numbered domains (defaults to
            90, matching the Hatteras Island hindcast).
        num_buffer_domains: Count of buffer domains padding each end
            (defaults to 15).
        first_gis_id: GIS ID of the first real domain (defaults to 1).
        domain_spacing_m: Alongshore width of one domain, in meters
            (defaults to 500).
    """

    num_real_domains: int = 90
    num_buffer_domains: int = 15
    first_gis_id: int = 1
    domain_spacing_m: float = 500.0

    @property
    def last_gis_id(self):
        """GIS ID of the last real domain."""
        return self.first_gis_id + self.num_real_domains - 1

    @property
    def total_domains(self):
        """Total padded array length (real domains + both buffers)."""
        return self.num_real_domains + 2 * self.num_buffer_domains

    @property
    def start_real_index(self):
        """Padded index of the first real domain (GIS first_gis_id)."""
        return self.num_buffer_domains

    @property
    def end_real_index(self):
        """Padded index one past the last real domain (exclusive)."""
        return self.num_buffer_domains + self.num_real_domains

    def gis_to_pad(self, gis_id):
        """Convert a GIS domain ID (e.g. 1-90) to its padded array index.

        Args:
            gis_id: GIS domain ID, or a numpy array/sequence of them.

        Returns:
            The padded array index (or array of indices).
        """
        return self.num_buffer_domains + (gis_id - self.first_gis_id)

    def pad_to_gis(self, pad_index):
        """Convert a padded array index back to its GIS domain ID.

        The inverse of `gis_to_pad`. Buffer domains have no GIS ID, so an
        index outside the real span returns a number below first_gis_id or
        above last_gis_id -- check the range yourself if it matters.

        Args:
            pad_index: Padded array index, or a numpy array/sequence of them.

        Returns:
            The GIS domain ID (or array of IDs).
        """
        return self.first_gis_id + (pad_index - self.num_buffer_domains)


DEFAULT_DOMAINS = DomainGeometry()
