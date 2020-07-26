# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the MIT License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import click

from XDOC.xdoc import xdoc
from XDOC import __version__


@click.command()
@click.option(
    "-i", "--i-otu", required=True, type=str,
    help="Input biom table."
)
@click.option(
    "-o", "--o-outdir", required=True, type=str,
    help="Output directory path."
)
@click.option(
    "-m", "--m-metadata", required=False, type=str,
    default=None, help="Metadata table."
)
@click.option(
    "-c", "--p-column", required=False, type=str,
    default=None, help="Column from metadata `-m` to use for "
                       "filtering based on values of `-v`."
)
@click.option(
    "-v", "--p-column-value", required=False, type=str, multiple=True,
    default=None, help="Filtering value to select samples based"
                       " on column passed to `-c`."
)
@click.option(
    "-q", "--p-column-quant", required=False, type=int,
    default=0, help="Filtering quantile / percentile for samples based on"
                    " column passed to `-c` (must be between 0 and 100)."
)
@click.option(
    "-fp", "--p-filter-prevalence", required=False, type=float,
    default=0, help="Filter features based on their minimum sample prevalence "
                    "(number >1 for sample counts: <1 for samples fraction)."
)
@click.option(
    "-fa", "--p-filter-abundance", required=False, type=float,
    default=0, help="Filter features based on their minimum sample abundance "
                    "(number >1 for abundance counts: <1 for abundance fraction)."
)
@click.option(
    "-f", "--p-filter-order", required=False, default='meta-filter',
    type=click.Choice(['meta-filter', 'filter-meta']),
    show_default=True, help="Order to apply the filters: 'filter-meta' first the prevalence/"
                            "abundance and then based on variable; 'meta-filter' first based "
                            "on variable and then the prevalence/abundance on the remaining."
)
@click.option(
    "-r", "--p-r", required=False, default=100, type=int,
    show_default=True, help="Number of bootstraps."
)
@click.option(
    "-subr", "--p-subr", required=False, default=0, type=int, show_default=True,
    help="If ignored will do bootstrap, alternatively an integer denoting size of subsample."
)
@click.option(
    "-pair", "--p-pair", required=False, default=None, type=str, multiple=True,
    help="One of the two names pair matching samples names."
)
@click.option(
    "-mov_avg", "--p-mov-avg", required=False, default=5, type=int, show_default=True,
    help="Moving average window to use for estimating where negative slope starts."
)
@click.option(
    "-ci", "--p-ci", required=False, type=float, default=(0.025, 0.5, 0.975,),
    multiple=True, show_default=True, help="Quantiles for confidence intervals."
)
@click.option(
    "-span", "--p-span", required=False, default=0.2, type=float,
    show_default=True, help="Span of loess smoothing."
)
@click.option(
    "-degree", "--p-degree", required=False, default=1., type=float, show_default=True,
    help="Degree of loess smoothing (If 1 linear, >1 polynomial)."
)
@click.option(
    "-family", "--p-family", required=False, show_default=True,
    default='symmetric', type=click.Choice(['gaussian', 'symmetric']),
    help="'gaussian' for least-squares fitting, 'symmetric' for robust fitting."
)
@click.option(
    "-iterations", "--p-iterations", required=False, default=4, type=int,
    show_default=True, help="Number of iterations for robust fitting."
)
@click.option(
    "-surface", "--p-surface", required=False, show_default=True,
    default='interpolate', type=click.Choice(['direct', 'interpolate']),
    help="'direct' estimation (slow exact) or 'interpolate' estimation (fast approximate)."
)
@click.option(
    "-cpus", "--p-cpus", required=False, default=1, type=int,
    show_default=True, help="Number of cores to use."
)
@click.option(
    "-nulls", "--p-nulls", required=False, default=1, type=int,
    show_default=True, help="Number of null models."
)
@click.option(
    "--non-zero/--no-non-zero", default=True,
    show_default=True, help="Only shuffle non zero values for null model."
)
@click.option(
    "--null/--no-null", default=False,
    show_default=True, help="Perform DOC on null models."
)
@click.option(
    "--verbose/--no-verbose", default=False
)
@click.version_option(__version__, prog_name="XDOC")
def standalone_xdoc(
        i_otu,
        o_outdir,
        m_metadata,
        p_column,
        p_column_value,
        p_column_quant,
        p_filter_prevalence,
        p_filter_abundance,
        p_filter_order,
        p_r,
        p_subr,
        p_pair,
        p_mov_avg,
        p_ci,
        p_span,
        p_degree,
        p_family,
        p_iterations,
        p_surface,
        p_cpus,
        p_nulls,
        non_zero,
        null,
        verbose
):

    xdoc(
        i_otu,
        o_outdir,
        m_metadata,
        p_column,
        p_column_value,
        p_column_quant,
        p_filter_prevalence,
        p_filter_abundance,
        p_filter_order,
        p_r,
        p_subr,
        p_pair,
        p_mov_avg,
        p_ci,
        p_span,
        p_degree,
        p_family,
        p_iterations,
        p_surface,
        p_cpus,
        p_nulls,
        non_zero,
        null,
        verbose
    )


if __name__ == "__main__":
    standalone_xdoc()
