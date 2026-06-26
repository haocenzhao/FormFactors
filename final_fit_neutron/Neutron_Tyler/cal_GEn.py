# Purpose:
#   - Evaluate GMn/(mu_n*GD) and its uncertainty at selected Q^2 points using Tyler's GetFF.
#   - Read user-provided mu_n*GE/GM (GEGM) and its uncertainty at the same Q^2 points.
#   - Compute GE/GD (GEGD) and its uncertainty by product error propagation (uncorrelated inputs).
#   - Enforce fixed numeric formats:
#       * Q2: 1 decimal
#       * GEGM, GEGM_Err: 3 decimals
#       * GMGD, GMGD_Err: computed then rounded to 3 decimals for downstream use
#       * GEGD, GEGD_Err: 3 decimals, computed using the rounded GMGD/GMGD_Err
#   - Write the output file with no trailing blank line.

from math import sqrt
import numpy as np

Q2_LIST = [2.9, 4.49, 4.51, 6.6, 9.9]

GEGM_INPUT = {2.9: 0.440, 4.49: 0.556, 4.51: 0.556, 6.6: 0.679, 9.9: 0.839}
GEGM_ERR_INPUT = {2.9: 0.028, 4.49: 0.078, 4.51: 0.048, 6.6: 0.086, 9.9: 0.143}


def GD(Q2: float) -> float:
    """Standard dipole form factor GD(Q2) = 1/(1+Q2/0.71)^2."""
    return 1.0 / (1.0 + Q2 / 0.71) ** 2


# Tyler's parameterization: GetFF(kID, Q2).
# For GMn: kID=4 returns (GMn/(mu_n*GD), d[GMn/(mu_n*GD)]).
from FF_newGMn_SBSproj import GetFF  # noqa: E402


def fmt_row(Q2, gd, GEGM, GEGM_Err, GMGD, GMGD_Err, GEGD, GEGD_Err) -> str:
    """Format one output line with fixed-width columns."""
    return (
        f"{Q2:6.1f}  "
        f"{gd:14.9e}  "
        f"{GEGM:7.3f}  {GEGM_Err:9.3f}  "
        f"{GMGD:7.3f}  {GMGD_Err:9.3f}  "
        f"{GEGD:7.3f}  {GEGD_Err:9.3f}"
    )


def main() -> None:
    out_dat = "GEnGMn_SBS.dat"

    rows = []
    for Q2 in Q2_LIST:
        gd = GD(Q2)

        GMGD_raw, GMGD_Err_raw = GetFF(4, Q2)
        GMGD = round(float(GMGD_raw), 3)
        GMGD_Err = round(float(GMGD_Err_raw), 3)

        GEGM = float(GEGM_INPUT[Q2])
        GEGM_Err = float(GEGM_ERR_INPUT[Q2])

        GEGD = round(GEGM * GMGD, 3)
        GEGD_Err = round(sqrt((GMGD * GEGM_Err) ** 2 + (GEGM * GMGD_Err) ** 2), 3)

        rows.append((Q2, gd, GEGM, GEGM_Err, GMGD, GMGD_Err, GEGD, GEGD_Err))

    # Header line only (line 1), then data (line 2+). No trailing blank line.
    header = (
        "# "
        + f"{'Q2[GeV^2]':>2}  "
        + f"{'GD':>10}  "
        + f"{'GEGM':>7}  {'GEGM_Err':>9}  "
        + f"{'GMGD':>7}  {'GMGD_Err':>9}  "
        + f"{'GEGD':>7}  {'GEGD_Err':>9}"
    )

    lines = [header]
    for r in rows:
        lines.append(fmt_row(*r))

    with open(out_dat, "w") as f:
        f.write("\n".join(lines))  # no extra newline at EOF

    # Print the same aligned table (also no extra blank line).
    print(f"[INFO] Wrote: {out_dat}")
    print("\n".join(lines))


if __name__ == "__main__":
    main()