# CODEX

## Purpose

This document describes how Codex should work in this project.

It is also intended as a compact example for other projects that want to be
easy for Codex, or another coding agent, to pick up safely.

## Read Order

Before making code changes, read:

1. `docs/DESIGN.md`
2. `docs/PHYSICS.md`
3. `docs/DECISIONS.md`
4. `docs/STATUS.md`
5. `docs/NEXT_TASK.md`

For implementation details, also read `docs/ALGORITHM.md`.

Use `docs/REFERENCES.md` when adding citations, validation context, or external
interfaces.

## Document Roles

- `DESIGN.md`: project identity, scope, and architecture
- `USAGE.md`: intended user workflows and package shape
- `PHYSICS.md`: authoritative physical model
- `ALGORITHM.md`: translation from physics to computation
- `DECISIONS.md`: accepted choices and deferred questions
- `STATUS.md`: current state and near-term milestone
- `ROADMAP.md`: staged future work
- `NEXT_TASK.md`: one concrete next task
- `REFERENCES.md`: bibliography and external tools
- `CODEX.md`: agent workflow and documentation conventions

## Working Rules

- Do not write code before checking the relevant design and physics documents.
- Keep the core generator-independent.
- Do not add a Pythia, HepMC, ROOT, or other external dependency to the core.
- Treat uncertain physics or event-record policy as an open question, not as an
  implicit implementation choice.
- Prefer small, reviewable changes that match the current milestone.
- Update documentation when a decision, scope boundary, or milestone changes.
- Keep `NEXT_TASK.md` focused on the immediate next task only.

## When Implementing

1. Confirm the requested change matches `STATUS.md` and `NEXT_TASK.md`.
2. Check `DECISIONS.md` for settled and deferred choices.
3. If the change affects physics, update or consult `PHYSICS.md` first.
4. If the change affects numerical behavior, update or consult `ALGORITHM.md`.
5. Implement the smallest coherent step.
6. Add or update tests for the behavior being changed.
7. Update `STATUS.md` when the project state changes.

## Handling Open Questions

If a task depends on an unsettled choice:

- add the question to `DECISIONS.md` or the relevant document,
- explain the options briefly,
- avoid silently choosing a physics or architecture policy unless the user has
  asked for that decision.

## Example Pattern For Other Projects

A Codex-friendly project should usually have:

- one stable design document,
- one domain-model document, if the project has specialized semantics,
- one algorithm or implementation-spec document,
- one decision log,
- one status file,
- one roadmap,
- one immediate next-task file,
- one agent guide like this file.

The important habit is separation of concerns: durable context belongs in stable
documents, current progress belongs in `STATUS.md`, and the next coding action
belongs in `NEXT_TASK.md`.
