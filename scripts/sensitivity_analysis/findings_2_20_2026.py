# =============================================================================
# WAVE CLIMATE SENSITIVITY ANALYSIS — HATTERAS ISLAND (1978–1997)
# SLR = 3.0 mm/yr | FLIP_SIGN = True
# =============================================================================
#
# RECOMMENDED CALIBRATED VALUES:
#
#   wave_period            = 8     # seconds  | LOW sensitivity — lines nearly
#                                  #             identical across 6–10s range;
#                                  #             fix to physically justified default
#
#   wave_asymmetry         = 0.7   # (0–1)    | LOW sensitivity — 0.6/0.7/0.8
#                                  #             indistinguishable island-wide;
#                                  #             fix to mid-range default
#
#   wave_height            = 2.0   # meters   | MODERATE sensitivity — spread
#                                  #             visible in mid-island erosion zone
#                                  #             (domains ~38–55) and southern section
#                                  #             (domains 75–90); 1.5m underpredicts,
#                                  #             2.5m overpredicts at south; 2.0m best
#                                  #             tracks observations and consistent
#                                  #             with WIS mean Hs offshore Hatteras
#
#   wave_angle_high_frac   = 0.2   # (0–1)    | MODERATE sensitivity — most
#                                  #             differentiation at northern boundary
#                                  #             (domains 1–5) and southern erosion
#                                  #             zone (domains 72–90); 0.4 underpredicts
#                                  #             southern erosion; 0.1 causes
#                                  #             over-accretion at north; 0.2 best
#                                  #             balances both ends
#
# NOTE: Wave climate parameters do NOT resolve local hotspot mismatches at
#   domains 63–67 (accretion peak) or 77–83 (Rodanthe erosion trough).
#   Those are driven by Oregon Inlet boundary conditions / sediment shadow
#   effects — a separate calibration problem, not wave climate.
#
# =============================================================================