"""Compatibility entry point for the superseded exploratory T19 script.

The submission-facing implementation now lives in
``src.baselines.run_multiview_fusion_benchmark`` and requires explicit input
paths so an external cohort cannot silently participate in training.
"""

from src.baselines.run_multiview_fusion_benchmark import main


if __name__ == "__main__":
    main()
