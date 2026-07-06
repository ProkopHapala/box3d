
# Box3D Visualization System Summary

## Overview

Box3D is a 3D rigid body physics engine (C17, by Erin Catto) with a sophisticated **sample testbed application** built on top of the **Sokol** header-only graphics framework. The testbed is a fully-featured 3D PBR renderer with shadows, ambient occlusion, image-based lighting, tonemapping, and debug overlays — all in pure C, cross-platform (OpenGL 4.5 / D3D11 / Metal).

CODEMAP:
Box3D Visualization System: Sokol-based 3D PBR Renderer
https://windsurf.com/codemaps/52148250-11bf-4d7f-92aa-ac962197f051-fe86ab10a43f3d18

## Sokol Framework

The externals live in [/home/prokophapala/git_SW/box3d/extern/sokol/](cci:9://file:///home/prokophapala/git_SW/box3d/extern/sokol:0:0-0:0):

- **[sokol_app.h](cci:7://file:///home/prokophapala/git_SW/box3d/extern/sokol/sokol_app.h:0:0-0:0)** (~580KB) — cross-platform window/event loop. Provides `sapp_desc`, [OnInit](cci:1://file:///home/prokophapala/git_SW/box3d/samples/main.cpp:75:0-117:1)/[OnFrame](cci:1://file:///home/prokophapala/git_SW/box3d/samples/main.cpp:331:0-432:1)/[OnEvent](cci:1://file:///home/prokophapala/git_SW/box3d/samples/main.cpp:119:0-309:1)/[OnCleanup](cci:1://file:///home/prokophapala/git_SW/box3d/samples/main.cpp:434:0-452:1) callbacks. The app entry point is [sokol_main()](cci:1://file:///home/prokophapala/git_SW/box3d/samples/main.cpp:473:0-531:1) at `@/home/prokophapala/git_SW/box3d/samples/main.cpp:474`.
- **[sokol_gfx.h](cci:7://file:///home/prokophapala/git_SW/box3d/extern/sokol/sokol_gfx.h:0:0-0:0)** (~1.2MB) — cross-backend GPU abstraction. Pipelines, buffers, images, passes, views, samplers — all via `sg_*` API.
- **[sokol_imgui.h](cci:7://file:///home/prokophapala/git_SW/box3d/extern/sokol/sokol_imgui.h:0:0-0:0)** — Dear ImGui integration.
- **[sokol_glue.h](cci:7://file:///home/prokophapala/git_SW/box3d/extern/sokol/sokol_glue.h:0:0-0:0)** — bridges `sokol_app` environment/swapchain to `sokol_gfx`.
- **[sokol_log.h](cci:7://file:///home/prokophapala/git_SW/box3d/extern/sokol/sokol_log.h:0:0-0:0)** — logging.

The implementation is compiled in a single TU at `@/home/prokophapala/git_SW/box3d/samples/gfx/sokol_impl.c:9-10`:
```c
#define SOKOL_IMPL
#include "sokol_gfx.h"
```

## Architecture Layers

### 1. Host Layer ([samples/host/](cci:9://file:///home/prokophapala/git_SW/box3d/samples/host:0:0-0:0))

- **[Camera](cci:2://file:///home/prokophapala/git_SW/box3d/samples/host/camera.h:42:0-255:1)** (`@/home/prokophapala/git_SW/box3d/samples/host/camera.h:43`) — orbit + fly camera. Dual-mode: Alt+drag orbits around a pivot, right-drag does FPS look with WASD. Produces `view`/`viewInv`/`proj`/`projInv` matrices together (no runtime inversion). Supports simulation→display transform (length unit scaling + Z-up reorientation). Double-precision eye (`m_worldEye`) for large-world rendering. Pick ray generation via `BuildPickRay()`.
- **`gui.cpp`** / `@/home/prokophapala/git_SW/box3d/samples/host/gui.h:38` — ImGui shell wrapping [sokol_imgui.h](cci:7://file:///home/prokophapala/git_SW/box3d/extern/sokol/sokol_imgui.h:0:0-0:0). `InitUI()` / `StartUIFrame()` / `RenderUI()` lifecycle. Edge-anchored panels via `BeginPanel()`.

### 2. Renderer ([samples/gfx/renderer.c](cci:7://file:///home/prokophapala/git_SW/box3d/samples/gfx/renderer.c:0:0-0:0) + [renderer.h](cci:7://file:///home/prokophapala/git_SW/box3d/samples/gfx/renderer.h:0:0-0:0))

The central orchestrator. `@/home/prokophapala/git_SW/box3d/samples/gfx/renderer.h:73` [InitRenderer()](cci:1://file:///home/prokophapala/git_SW/box3d/samples/gfx/renderer.c:645:0-1273:1) creates all GPU resources. The frame pipeline (see `@/home/prokophapala/git_SW/box3d/samples/gfx/renderer.c:2594`) is:

1. **[PreSceneWork](cci:1://file:///home/prokophapala/git_SW/box3d/samples/gfx/renderer.c:2102:0-2127:1)** (`renderer.c:2104`):
   - Upload all per-frame instance data (cubes, spheres, capsules, meshes)
   - Rebuild IBL if sun/turbidity changed (`RebuildImageBasedLightingIfDirty`)
   - Render 3-cascade shadow maps ([RenderShadowCascades](cci:1://file:///home/prokophapala/git_SW/box3d/samples/gfx/renderer.c:2073:0-2100:1))
   - **Depth pre-pass** — rasterize analytic shape depth into R32F linear-depth target
   - **GTAO** — `PrefilterDepth` → `ComputeNoisyResult` → `Denoise` (XeGTAO port from Esoterica)

2. **[DrawSceneIntoHdr](cci:1://file:///home/prokophapala/git_SW/box3d/samples/gfx/renderer.c:2129:0-2152:1)** (`renderer.c:2134`) — 4x MSAA RGBA16F scene target:
   - Opaque cubes, spheres (impostor shader), capsules (impostor shader), triangle geometry (hulls/meshes/heightfields)
   - Procedural Preetham sky (fullscreen triangle at far plane)
   - Edge overlay (convex/concave/flat edge classification, MSAA)

3. **[DrawTransparentIntoResolve](cci:1://file:///home/prokophapala/git_SW/box3d/samples/gfx/renderer.c:2192:0-2536:1)** (`renderer.c:2194`) — back-to-front sorted transparent shapes blended over the resolved HDR scene

4. **Transparent edge overlay** — edges on top of transparent blend

5. **Highlight mask** — R8 mask rasterized for hovered/selected shapes

6. **Present pass** (swapchain):
   - **AgX tonemap** — HDR → display sRGB with exposure control
   - **Highlight outline** — silhouette compositing from the mask
   - **Overlay lines/points** — post-tonemap (display-referred, WYSIWYG colors), with occlusion modes (HIDE/DIM/DASHED) sampling the GTAO linear-depth

Key renderer types:
- `FrameInput` (`@/home/prokophapala/git_SW/box3d/samples/gfx/renderer.h:19`) — per-frame camera matrices, draw origin, debug mode, shadow/AO toggles, Z-up flag
- `Sun` (`renderer.h:65`) — directional light with color, ambient, strength
- `RenderStats` (`renderer.h:107`) — per-frame draw call / instance counts
- [RenderFrame()](cci:1://file:///home/prokophapala/git_SW/box3d/samples/gfx/renderer.c:2593:0-2673:1) / [RenderFrameOffscreen()](cci:1://file:///home/prokophapala/git_SW/box3d/samples/gfx/renderer.c:2675:0-2735:1) — two entry points (swapchain vs. caller-owned target)

### 3. Rendering Subsystems ([samples/gfx/](cci:9://file:///home/prokophapala/git_SW/box3d/samples/gfx:0:0-0:0))

Each is a self-contained C module with `Init`/`Shutdown` and per-frame submit:

- **`scene_target`** (`@/home/prokophapala/git_SW/box3d/samples/gfx/scene_target.h:40`) — 4x MSAA RGBA16F color + D32F depth + single-sample resolve. `SCENE_SAMPLE_COUNT = 4`.
- **`shadow`** (`@/home/prokophapala/git_SW/box3d/samples/gfx/shadow.h:45`) — 3-cascade PSSM, 2048² D32F array, texel-snapped frustum fitting, PCF linear sampler.
- **`gtao`** (`@/home/prokophapala/git_SW/box3d/samples/gfx/gtao.h:54`) — XeGTAO port. Depth pre-pass → 5-mip prefiltered depth chain → main trace (9 slices × 3 steps at Ultra) → bilateral denoise (up to 5 passes). Quality tiers: Medium/High/Ultra.
- **`ibl`** (`@/home/prokophapala/git_SW/box3d/samples/gfx/ibl.h:43`) — Karis split-sum BRDF LUT (256²), prefiltered sky cubemap (256², 7 mips), 9-band spherical harmonics for diffuse. Rebuilds on sun/turbidity change.
- **`sky`** (`@/home/prokophapala/git_SW/box3d/samples/gfx/sky.h:27`) — Preetham analytic daylight, fullscreen triangle at reverse-Z far plane, below-horizon fade.
- **`tone_map`** (`@/home/prokophapala/git_SW/box3d/samples/gfx/tone_map.h:29`) — AgX tonemap with exposure stops and saturation control.
- **`overlay`** (`@/home/prokophapala/git_SW/box3d/samples/gfx/overlay.h:35`) — SDF-coverage-AA lines and points, post-tonemap, depth-tested via GTAO linear-depth. Occlusion modes: HIDE, DIM (reduced alpha), DASHED.
- **`edges`** (`@/home/prokophapala/git_SW/box3d/samples/gfx/edges.h:34`) — Per-shape edge rendering for hulls/meshes/heightfields. Convex/concave/flat classification. Opaque edges in MSAA pass, transparent edges in separate single-sample pass.
- **`highlight_mask`** + **`highlight_outline`** — Two-pass selection highlight: R8 mask rasterization → fullscreen outline compositing. Hover (0.5) vs. select (1.0) states.
- **`geometry_registry`** (`@/home/prokophapala/git_SW/box3d/samples/gfx/geometry_registry.h:52`) — Ref-counted GPU buffer store for triangle meshes. Keyed by content hash. `FindMesh`/`RegisterMesh`/`AddMeshReference`/`ReleaseMeshReference` lifecycle. Per-frame instance arenas, one shared GPU storage buffer, per-draw `inst_base` uniform offset.
- **`debug_shapes`** (`@/home/prokophapala/git_SW/box3d/samples/gfx/debug_shapes.h:35`) — Bridge from Box3D collision shapes (`b3HullData`, `b3MeshData`, `b3HeightFieldData`) to the geometry registry. Flat-shaded (duplicated vertices per face).
- **`debug_adapter`** (`@/home/prokophapala/git_SW/box3d/samples/gfx/debug_adapter.h:19`) — Box3D→renderer bridge. `AttachToWorldDef()` hooks debug draw callbacks. `MakeDebugDraw()` returns a `b3DebugDraw` struct. Per-shape material override, ground grid tagging, compound child culling, selection state.
- **`text`** (`@/home/prokophapala/git_SW/box3d/samples/gfx/text.h:37`) — World-space and screen-space text labels. Accumulated per-frame, drained by the ImGui shell via `GetTextAt()`.
- **`utility`** (`@/home/prokophapala/git_SW/box3d/samples/gfx/utility.h:25`) — `Vec4`, `Mat4` (column-major), [MakePerspective](cci:1://file:///home/prokophapala/git_SW/box3d/samples/gfx/utility.h:85:0-103:1) (reverse-Z, 0..1 clip), [MakeLookAt](cci:1://file:///home/prokophapala/git_SW/box3d/samples/gfx/utility.h:150:0-163:1), [MakeViewAndInverse](cci:1://file:///home/prokophapala/git_SW/box3d/samples/gfx/utility.h:165:0-189:1) (no runtime inversion), [MakeMat4FromTransform](cci:1://file:///home/prokophapala/git_SW/box3d/samples/gfx/utility.h:191:0-202:1), sRGB→linear color conversion.

### 4. Draw API ([samples/gfx/draw.h](cci:7://file:///home/prokophapala/git_SW/box3d/samples/gfx/draw.h:0:0-0:0))

The high-level drawing interface used by samples (`@/home/prokophapala/git_SW/box3d/samples/gfx/draw.h:27`):

- **Shapes**: `DrawCube`, `DrawSphere`, `DrawCapsule` (+ `Ex` variants with PBR material params), `DrawSolidSphere`, `DrawSolidCapsule`, `DrawHull`, `DrawPlane`
- **Overlays**: `DrawLine`, `DrawPoint`, `DrawArrow`, `DrawCross`, `DrawAabb`, [DrawBounds](cci:1://file:///home/prokophapala/git_SW/box3d/samples/host/camera.h:171:1-179:2), `DrawAxes`, `DrawGrid`, `DrawGroundGrid`, `DrawTriangle`, `DrawWireSphere`, `DrawWireCapsule`
- **Text**: `DrawString3D` (printf-style, world-space)
- **Origin**: `SetDrawOrigin` / `GetDrawOrigin` — double-precision origin shift so float coordinates stay small far from world origin

All calls demote against the draw origin in double precision before handing float coordinates to the renderer.

### 5. Sample Framework ([samples/sample.h](cci:7://file:///home/prokophapala/git_SW/box3d/samples/sample.h:0:0-0:0), [sample.cpp](cci:7://file:///home/prokophapala/git_SW/box3d/samples/sample.cpp:0:0-0:0))

`@/home/prokophapala/git_SW/box3d/samples/sample.h:107` — [Sample](cci:2://file:///home/prokophapala/git_SW/box3d/samples/sample.h:106:0-229:1) base class:

- Virtual [Step()](cci:1://file:///home/prokophapala/git_SW/box3d/samples/sample_shapes.cpp:794:1-822:2) (physics), [Render()](cci:1://file:///home/prokophapala/git_SW/box3d/samples/sample.h:121:1-121:25) (drawing), [DrawControls()](cci:1://file:///home/prokophapala/git_SW/box3d/samples/sample_shapes.cpp:773:1-792:2) (ImGui panel), [Keyboard()](cci:1://file:///home/prokophapala/git_SW/box3d/samples/sample.h:173:1-175:2), `MouseDown/Up/Move()`
- `CreateWorld()` sets up the Box3D world with debug draw callbacks
- `AddGroundBox()` — static ground with procedural grid material
- Recording/replay support (`StartRecording`/`FinishRecording`)
- [SampleContext](cci:2://file:///home/prokophapala/git_SW/box3d/samples/sample.h:22:0-104:1) (`sample.h:23`) — shared state: camera, window size, simulation settings (hertz, substeps, worker count), render toggles (shadows, GTAO, Z-up), UI state
- [RegisterSample()](cci:1://file:///home/prokophapala/git_SW/box3d/samples/sample.h:248:0-248:83) — static registration pattern, ~30+ samples across categories (Shapes, Collision, Joints, Character, etc.)

### 6. Main Loop ([samples/main.cpp](cci:7://file:///home/prokophapala/git_SW/box3d/samples/main.cpp:0:0-0:0))

`@/home/prokophapala/git_SW/box3d/samples/main.cpp:332` [OnFrame()](cci:1://file:///home/prokophapala/git_SW/box3d/samples/main.cpp:331:0-432:1):

1. Camera update (dt, window size)
2. `SetDrawOrigin(camera.DrawOrigin())` — sync render origin
3. [ResetFrameArena()](cci:1://file:///home/prokophapala/git_SW/box3d/samples/gfx/renderer.c:2928:0-2952:1) — clear all per-instance arenas
4. [sample->Step()](cci:1://file:///home/prokophapala/git_SW/box3d/samples/sample_shapes.cpp:794:1-822:2) — physics + HUD text
5. [sample->Render()](cci:1://file:///home/prokophapala/git_SW/box3d/samples/sample.h:121:1-121:25) — fill instance/overlay arenas via `Draw*` calls and `b3DebugDraw`
6. Build `FrameInput` from camera state
7. [RenderFrame(&sc, &fi)](cci:1://file:///home/prokophapala/git_SW/box3d/samples/gfx/renderer.c:2593:0-2673:1) — full GPU pipeline
8. `StartUIFrame(dt)` — ImGui + text drain
9. `RenderUI(&sc)` — ImGui pass on swapchain
10. `sg_commit()`
11. Software 60 Hz frame limiter ([LimitFrameRate](cci:1://file:///home/prokophapala/git_SW/box3d/samples/main.cpp:311:0-329:1))

## Key Design Patterns

- **Header-only Sokol** — single implementation TU, backend selected by CMake target
- **C renderer, C++ samples** — renderer is pure C (`extern "C"`), samples are C++ inheriting from [Sample](cci:2://file:///home/prokophapala/git_SW/box3d/samples/sample.h:106:0-229:1)
- **Per-frame arena** — all instance data accumulates in CPU arenas, uploaded once per frame, reset at frame start
- **No runtime matrix inversion** — view/proj and their inverses are built together analytically
- **Reverse-Z** — near/far swapped in perspective, GL 4.5 required for `glClipControl`
- **Double-precision origin** — `b3Pos` (double3) for draw origin, demoted to float only at the renderer boundary
- **Host-agnostic** — [RenderFrameOffscreen()](cci:1://file:///home/prokophapala/git_SW/box3d/samples/gfx/renderer.c:2675:0-2735:1) allows non-sokol_app hosts (e.g. GLFW)
- **Shaders** — compiled via `sokol-shdc` (see [shaders/](cci:9://file:///home/prokophapala/git_SW/box3d/samples/shaders:0:0-0:0) directory), generated headers included in [renderer.c](cci:7://file:///home/prokophapala/git_SW/box3d/samples/gfx/renderer.c:0:0-0:0)

## Using This as a Template

The renderer ([samples/gfx/](cci:9://file:///home/prokophapala/git_SW/box3d/samples/gfx:0:0-0:0)) is cleanly separated from the physics engine and could be extracted as a standalone 3D debug renderer. The key dependencies are:

- **Sokol** headers (gfx, app, glue, imgui, log)
- **Box3D** math types (`b3Vec3`, `b3Transform`, `b3Quat`, etc.) — but these could be swapped
- The [utility.h](cci:7://file:///home/prokophapala/git_SW/box3d/samples/gfx/utility.h:0:0-0:0) `Mat4`/`Vec4` types are self-contained and designed to be cut-pasted

The architecture is modular: each rendering subsystem (shadow, GTAO, IBL, sky, tonemap, overlay, edges, highlights) has its own `Init`/`Shutdown` and per-frame submit function, making it easy to add, remove, or replace individual features.