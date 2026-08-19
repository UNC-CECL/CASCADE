"""Identifying info for one completed CASCADE run.

Bundles the handful of values every plotting function in cascade_pipeline needs
(run name/dir, period, wave height, sign convention) so call sites pass one
object instead of five loose scalars pulled from module globals.
"""

import dataclasses


@dataclasses.dataclass(frozen=True)
class RunInfo:
    """Identifying info for one completed CASCADE run.

    Attributes:
        run_name: Run name used in output filenames (e.g. run_name_hs).
        run_dir: Output directory for this run's files.
        start_year: Calendar year the run starts (period START_YEAR).
        end_year: Calendar year the run ends.
        Hs: Fixed significant wave height used for this run (m), shown in
            figure titles/legends. None omits it.
        flip_sign_model: CASCADE's x_s_TS increases landward (erosion).
            True (the convention used throughout cascade_pipeline's plotting)
            flips model output so a larger value means more seaward.
        background_erosion_on: Whether DOMAIN_BE_RATES has any non-zero
            entries for this run; used only for the "BE=on/off" label in
            figure titles and the GIF caption.
        model_name: Model name shown in figure titles/captions, e.g.
            "CASCADE". cascade_pipeline.shoreline is itself CASCADE-specific
            (it reads cascade.barrier3d / b3d.x_s_TS directly), so this
            defaults to "CASCADE" rather than a fully generic placeholder --
            override it if you're comparing against a different model run
            through the same figures.
    """

    run_name: str
    run_dir: str
    start_year: int
    end_year: int
    Hs: float = None
    flip_sign_model: bool = True
    background_erosion_on: bool = True
    model_name: str = "CASCADE"
