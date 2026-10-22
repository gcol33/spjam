# Story: spjam

## Claim
Joint spatial models that fit the ecological process and the sampling process together, coupled through a shared latent field, so the ecological map is estimated after preferential sampling is accounted for.

## Audience
Ecologists building species maps from opportunistic or effort-biased data, familiar with spatial GLMMs or GAMs.

## Relation
- to: none (positions itself against standard spatial GLMMs/GAMs and, by analogy, occupancy models)
- phrase: none
- home: README (Standard spatial model or spjam?)

## Arcs

### A1 Sampling as its own process
- question: does my map show the species or where people looked?
- answer: the sampling intensity is a second spatial process sharing a field with the ecology; beta_share measures the coupling
- evidence: ?spjam, tests
- home: README (Modelling where observers sampled); no article yet
- also: DESCRIPTION

### A2 Shared versus independent fit
- answer: latent_shared() against latent_independent(); the difference is the correction
- home: README

### A3 Building blocks
- answer: families, spatial fields, coupling, sampling models, priors
- home: help pages

## Non-goals
None stated.

## Documents
- DESCRIPTION: claim
- README: claim, supports list, install (compiled), first fit, coupling section, comparison table
- articles: none; no pkgdown site

## Open
- No _pkgdown.yml and no site; README links no reference pages and has no .bib/.ris download line.
- HMC is the spjam() default backend while the README calls Laplace the supported path.
