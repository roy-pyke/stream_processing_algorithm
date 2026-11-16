# Reference material (not redistributed)

The papers and third-party code we worked from are copyrighted and are not part of this repository. This
list says where to get each of them. The full bibliography is at the end of `paper/paper.tex`.

| Local file (ignored, `refs/`) | Work | Where to get it |
|---|---|---|
| `omnisketch_arxiv.pdf` | W. R. Punter, O. Papapetrou, M. Garofalakis. OmniSketch: Efficient multi-dimensional high-velocity stream analytics with arbitrary predicates. PVLDB 17(3), 2023 | <https://arxiv.org/abs/2309.06051> |
| `omnisketch_sigrec.pdf` | Punter, Papapetrou, Garofalakis. OmniSketch: Streaming data analytics with arbitrary predicates. SIGMOD Record 54(1), 2025 | SIGMOD Record, ACM Digital Library |
| `omnisketch_vldbj.pdf` (missing) | Punter, Papapetrou, Garofalakis. OmniSketch: Multi-dimensional update stream analytics with arbitrary predicates. The VLDB Journal 35(1), 2026 | open access, <https://doi.org/10.1007/s00778-025-00960-6> (the publisher blocks scripted downloads; use a browser) |
| `dbdbd25.pdf` | Punter, Papapetrou, Garofalakis. Extended abstract of the journal version, DBDBD 2025 | DBDBD 2025 proceedings |
| `vldbphd18.pdf` | W. R. Punter. Towards flexible self-tuning data stream management systems. VLDB PhD Workshop 2018 | VLDB PhD Workshop proceedings |
| `A Bird's-Eye View of Complex Streaming Data Analytics.pdf` | M. Garofalakis. A bird's-eye view of complex streaming data analytics. DEBS 2025 | <https://doi.org/10.1145/3701717.3733848> |
| `hydra.pdf` | A. Manousis et al. Enabling efficient and general subpopulation analytics in multidimensional data streams. PVLDB 15(11), 2022 | <https://arxiv.org/abs/2208.04927> |
| `cohen_multiobjective.pdf` | E. Cohen. Multi-objective weighted sampling. HotWeb 2015 | <https://arxiv.org/abs/1509.07445> |
| `ting_ats.pdf` | D. Ting. Adaptive threshold sampling. SIGMOD 2022 | <https://arxiv.org/abs/1708.04970> |
| `ting_uss.pdf` | D. Ting. Data sketches for disaggregated subset sum and frequent item estimation. SIGMOD 2018 | <https://arxiv.org/abs/1709.04048> |
| `svoila.pdf` | T. D. Nguyen et al. Stratified random sampling from streaming and stored data (S-VOILA). EDBT 2019 | <https://arxiv.org/abs/1801.09039> |
| `cvopt.pdf` | T. D. Nguyen et al. Random sampling for group-by queries (CVOPT). ICDE 2020 | <https://arxiv.org/abs/1909.02629> |
| `omni_src/*.java` | the authors' OmniSketch implementation | <https://github.com/wiegerpunter/omnisketch>, commit `8dc67b5`; fetched and built unmodified by `scripts/setup_upstream.sh` |
| (not stored) | the journal version's code, configuration files and analysis notebooks | <https://github.com/wiegerpunter/omni-deletes>, commit `17b969e`; fetched and built unmodified by `scripts/setup_omni_deletes.sh` (bounded deletions, #25) |
| (not stored) | the journal version's data pipeline (generator, residue/noise split, mixing of inserts and deletes) | <https://github.com/wiegerpunter/dataset-pipeline>, commit `06d5dc1`; re-implemented by `DeletesHarness gen` and `bin/mkdelstream` |
