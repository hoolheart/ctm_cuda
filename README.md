ctm_cuda
=============

Implementation of CTM by CUDA

### lane-intersection structure

The previous version builds the model directly on the cell-link level, which makes it hard to scale.
To solve it, the new version will focus on the lane-intersection level.
During the operations on the lanes and the intersections, the related cells and links are automatically created or modified.

Another difference with the previous version is that the `access` field is transferred from `cell` to `link`,
so that the model can handle more complicated phases.

The difference between this repository and the one entitled `hoolheart/CTM_C_Program.git` is that this repository runs the simulation on GPU instead of CPU.

