# rpq9MultiSoftMod Reference


## rpq9MultiSoftMod Node

### Overview
The rpq9MultiSoftMod node deforms one or more geometries, similar to a softMod connected in parallel.
The more softMod calculations you add, the shorter the processing time will be compared to a standard softMod.

### Attributes

| Long              | Short | Type       | Default  | Description                      |
| ----------------- | ----- | ---------- | --------| -------------------------------- |
| `inputData`       | `id`  | `compound` |          |The data from the index where this attribute exists will be used in the calculation. |
| `localEnvelope`   | `le`  | `float`    | 1.0      | Envelope for each softMod calculation. |
| `centerMatrix`    | `cm`  | `matrix`   | identity | Center matrix for each softMod calculation. |
| `modifyMatrix`    | `mm`  | `matrix`   | identity | Modify matrix for each softMod calculation. |
| `falloffRadius`   | `fr`  | `float`    | 1.0      | Deformation radius for each softMod calculation. |
| `localWeightList` | `lwl` | `compound` |          | Vertex weight lists for each softMod calculation. |
| `localWeights`    | `lw`  | `float`    | 1.0      | Vertex weights for each softMod calculation. |
| `falloffMode`     | `fm`  | `enum`     | 2        | falloff calculation mode. Please switch between 0 (None), 1 (Linear), 2 (Smooth), 3 (Spline), 4 (SmoothStep), and 5 (EaseInOut) depending on your purpose. |

---
<br>


## Commands


```python
import maya.cmds
import rpq9MultiSoftMod

# create rpq9MultiSoftMod to selection like standard deformers.
cmds.CreateRpq9MultiSoftMod

# show rpq9MultiSoftMod options like standard deformers.
cmds.CreateRpq9MultiSoftModOptions

# show window for setting deformer weights.
rpq9MultiSoftMod.deformerWeightUI.showWindow()

# For other commands, please check the files under ./scripts/rpq9MultiSoftMod.
# Frequently used functions are defined in "core.py"
```


## UI

### - deformerWeightUI
#### Overview
This is an auxiliary window for setting local weights in the same way as the standard deformer weight setup procedure.
as well as for reading and writing deformer weights that include local weights.

You can copy values between standard deformer weights and local weight values bidirectionally.

You can also read and write deformer weight data that includes local weights.

#### Instructions for local weight settings

1. Select the geometry with rpq9MultiSoftMod applied.

2. Press the `Set` button. If there are multiple rpq9MultiSoftMod, select the name of the deformer whose weights you want to set from the Deformer menu.

3. Press the `Paint Tool` button.

4. Edit deformer weights.

5. Set the index for the local weight to be configured in `localIndex`. (If the connections are simple, you can select the controller and use Set Sel to retrieve the index.)

6. Press the `->` button.

7. Press the `clear Weight` button.

#### Others

The `stash weight` button can be used to temporarily save the standard deformer weights.
The `unstash weight`”` button restores the weights that were previously stashed.

The `import` and `export` buttons allow you to read and write deformer weights for the currently specified deformer.