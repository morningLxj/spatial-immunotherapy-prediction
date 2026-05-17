# PLOS ONE Public `main` Realignment Design

## Goal

Realign the public GitHub `main` branch of `morningLxj/spatial-immunotherapy-prediction` to match the current PLOS ONE manuscript positioning.

The public repository should present a conservative, manuscript-consistent narrative focused on:

- data accession and processing transparency
- feature selection and genetic prioritization
- spatial pattern analysis
- computational consistency analysis
- reporting and reproducibility support

The repository should no longer read like a causal-claim or clinical-deployment package.

## Problem Summary

The current `main` branch still exposes an older public narrative:

- title foregrounds `causal inference`
- README emphasizes `final prognostic model`, `clinical and external validation`, and `MR-prioritized genes`
- repository layout is documented as `code/`, `data/`, `docs/`, `results/`
- manuscript-facing `Code Availability` language reportedly expects a public release layout centered on `01_data_accession_manifest` through `06_session_info`

This creates a mismatch between:

- the current PLOS ONE manuscript
- the cover-letter positioning
- the public code repository cited in the manuscript

## Constraints

- Do not rewrite analysis code unless needed for documentation accuracy.
- Do not physically rename the entire repository tree in this pass.
- Do not expose private data, raw patient-level files, or manuscript binaries.
- Do not preserve older stronger claims in public-facing text if they conflict with the current manuscript.
- Keep implementation low-risk and reversible.

## Approaches Considered

### Approach 1: Full Physical Restructure

Rename the repository into manuscript-style top-level directories such as `01_data_accession_manifest` through `06_session_info`.

**Pros**

- Matches manuscript wording exactly.
- Minimizes visible discrepancy for readers browsing the root.

**Cons**

- Highest implementation risk.
- Likely breaks links, script paths, and existing documentation.
- Requires broader regression checking than is justified at this stage.

### Approach 2: README-Only Refresh

Keep the repository structure as-is and update only the README plus a few docs.

**Pros**

- Fastest option.
- Lowest technical risk.

**Cons**

- Root directory still visibly conflicts with the manuscript's described public layout.
- Leaves ambiguity for editors or reviewers who inspect the repository.

### Approach 3: Public Release View Over Existing Structure

Keep the physical directories (`code/`, `data/`, `docs/`, `results/`) but introduce a manuscript-facing public release view that maps the current repository into the PLOS ONE structure and updates all public text to the conservative manuscript framing.

**Pros**

- Best balance of consistency and safety.
- Avoids breaking code paths.
- Gives the manuscript a stable, defensible `main` branch target.

**Cons**

- Requires careful wording so the mapping feels intentional, not improvised.
- Public readers see both the physical tree and the manuscript-facing mapped view.

**Recommendation**

Use Approach 3.

## Scope

### In Scope

- rewrite the root `README.md`
- create a manuscript-facing public release layout document
- add lightweight top-level documentation entry points corresponding to the manuscript's `01` to `06` structure
- align repository wording with the current PLOS ONE manuscript tone
- update code/data availability-facing wording inside repository docs where needed
- keep the repository easy to browse for editors and reviewers

### Out of Scope

- major code refactoring
- re-running analyses
- changing figure/table outputs
- full directory renaming
- adding new scientific claims
- exposing large or restricted datasets

## Target Narrative

The public `main` branch should communicate:

- this is a curated public release for the current manuscript
- the repository supports reproducibility and transparency, not clinical deployment
- MR is used as genetic prioritization support, not as proof of biological causality
- external validation is presented as supportive transportability assessment, not product-grade prediction
- computational attenuation analyses are internal consistency checks, not functional perturbation experiments

The public `main` branch should avoid foregrounding:

- causal-inference branding in the title
- strong clinical prediction framing
- deployable-model implications
- mechanistic overstatement

## Information Architecture

### Physical Structure

Keep the current repository tree intact:

- `code/`
- `data/`
- `docs/`
- `results/`

### Manuscript-Facing Public Release View

Introduce a documented mapping from the current tree into six public-release sections:

1. `01_data_accession_manifest`
2. `02_processing_and_qc`
3. `03_feature_selection_and_genetic_prioritization`
4. `04_spatial_and_consistency_analysis`
5. `05_reporting_assets`
6. `06_session_info`

Each section should point to existing files or directories rather than duplicating analysis assets.

## Planned Artifacts

### 1. Root README Rewrite

Rewrite `README.md` so it:

- uses a conservative PLOS ONE-compatible title
- leads with study transparency and analytical overview
- describes the repository as a curated public release
- links readers first to the six public-release sections
- de-emphasizes model superiority and clinical utility language
- removes older phrasing such as:
  - `causal inference-guided`
  - `final prognostic model`
  - `clinical and external validation`
  - `MR-prioritized genes` as a headline claim

### 2. Public Release Layout Document

Create a new document, likely `docs/public_release_layout.md`, that:

- explains why the public repository is presented through six manuscript-facing sections
- maps each section to actual files/directories in the repository
- helps editors and reviewers navigate without needing the private local workspace

### 3. Lightweight Top-Level Entry Documents

Add six lightweight top-level markdown files:

- `01_data_accession_manifest.md`
- `02_processing_and_qc.md`
- `03_feature_selection_and_genetic_prioritization.md`
- `04_spatial_and_consistency_analysis.md`
- `05_reporting_assets.md`
- `06_session_info.md`

Each file should be short and act as a stable entry point for manuscript readers.

These files should:

- summarize what the section contains
- link to the relevant directories/files
- avoid duplicating detailed content that already exists elsewhere

### 4. Supporting Docs Refresh

Adjust existing docs as needed so they no longer conflict with the new public-release framing:

- `docs/repository_layout.md`
- `docs/data_notes.md`
- `docs/latest_submission_sync.md`
- `docs/reproducibility_checklist.md`

The aim is alignment, not complete replacement.

### 5. Citation and Metadata Consistency

Review `CITATION.cff` and any repository-facing citation text to ensure the title and message do not overstate causal or clinical claims relative to the current manuscript.

## Content Rules

### Preferred Terms

- `genetic prioritization`
- `supportive evidence`
- `spatial pattern analysis`
- `computational consistency analysis`
- `external support`
- `transportability assessment`
- `curated public release`

### Terms To Minimize or Remove

- `causal inference-guided` in the public-facing title
- `causal anchor`
- `final prognostic model` as a central claim
- `clinical utility`
- `deployable prediction`
- `mechanistic` framing unless explicitly qualified

## Data Flow

The public reader flow should be:

1. open `README.md`
2. understand the conservative study framing
3. navigate via the six public-release entry points
4. drill into underlying `code/`, `docs/`, or `results/` references as needed

This means `README.md` and the six entry files become the public navigation layer, while the current repository tree remains the implementation layer.

## Error Handling and Risk Control

- Do not move or rename existing analysis directories in this pass.
- Do not edit scientific claims beyond aligning them downward to the current manuscript framing.
- Preserve working script paths.
- Avoid dead links by preferring relative links to already existing files.
- If any existing doc strongly conflicts with the new framing and cannot be safely reconciled, replace its wording rather than layering contradictory statements.

## Validation Plan

Before calling the update complete:

- verify the README title and overview reflect the conservative manuscript framing
- verify each of the six public-release entry files exists and links correctly
- verify the public release layout doc clearly maps manuscript sections to repository assets
- verify no public-facing top-level doc still foregrounds old causal/clinical-deployment language
- verify links resolve within GitHub-style relative navigation
- verify `git status` is clean apart from intentional files before final commit

## Implementation Sequence

1. rewrite `README.md`
2. add `docs/public_release_layout.md`
3. add the six top-level entry markdown files
4. refresh the most visible conflicting docs
5. review citation/metadata wording
6. run a final consistency pass across public-facing repository text

## Success Criteria

This redesign succeeds if:

- the public `main` branch reads as a PLOS ONE-consistent public release
- editors can understand the repository without seeing the private workspace
- manuscript `Code Availability` can credibly point to `main`
- the repository no longer reads like a stronger causal or clinical-prediction package than the manuscript itself
