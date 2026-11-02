:orphan:

.. _rda-case-details:

Acquisition tests: case details
===============================

Supplementary evidence for :doc:`remote_data_acquisition`, sections 3 and 4.

**Case labels.** The acquisition tests (D-Bench, the project's acquisition benchmark) used ten historical South
Atlantic cyclone trajectories of 1979–1980. Each trajectory is identified by an internal track label (for example
``19791150``). The labels exist only for traceability: they are not dates and are not parameters you pass to the
program. The same labels identify the cases in the run table
``_static/remote_data_acquisition/stageB_runs.csv`` (produced by ``benchmarks/d_bench/stageb_report.py``) and in
the permanent records under ``docs/evolution/``.

Energetics of NCAR files against the CDS
----------------------------------------

Summarised in :ref:`section 4.2 <rda-science>`. Tolerance per column: max(10\ :sup:`-4` × column scale,
10\ :sup:`-9`), fixed before the comparison. The same verdict applies to every NCAR method in a case, because the
fields are identical.

.. list-table:: Energetics of NCAR files against the CDS (37 levels)
   :header-rows: 1
   :widths: 16 18 20 46

   * - Case
     - Verdict
     - Largest difference / tolerance
     - Columns outside the tolerance
   * - 19790205
     - equivalent
     - 1.6 × 10\ :sup:`-5`
     - —
   * - 19790535
     - equivalent
     - 7.1 × 10\ :sup:`-5`
     - —
   * - 19790585
     - equivalent
     - 7.7 × 10\ :sup:`-5`
     - —
   * - 19800296
     - equivalent
     - 0.66
     - —
   * - 19791008
     - **not equivalent**
     - 1.24
     - ∂Az/∂t
   * - 19790644
     - **not equivalent**
     - 3.40
     - ∂Az/∂t, RGz
   * - 19800256
     - **not equivalent**
     - 3.48
     - BAz, ∂Az/∂t, ∂Ae/∂t, RGz, RGe
   * - 19790166
     - **not equivalent**
     - 5.18
     - Az, BAz, ∂Az/∂t, ∂Ae/∂t, RGz
   * - 19791150
     - **not equivalent**
     - 6.66
     - Az, Ae, Ca, BAz, Gz, ∂Az/∂t, ∂Ae/∂t, RGz, RGe
   * - 19800318
     - no CDS file
     - —
     - identical between the NCAR methods

- The temperature-substitution diagnosis was made on case ``19790644``.
- The largest temperature difference alone does not predict the verdict: ``19790205`` and ``19800296`` are
  equivalent with differences of the same order.
- The 13-level ARCO comparison and the vertical-resolution analysis (:ref:`rda-vertical`) used case ``19791150``.
