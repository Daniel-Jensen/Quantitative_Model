#!/usr/bin/env bash
# PreToolUse(Bash) gate: block `git commit` of model/code changes unless the
# living docs (STATE.md, PROCESS.md, HANDOFF.md) are updated in the SAME commit.
#
# Project policy: keep the docs in lockstep with the code — not just CLAUDE.md.
# Reads the tool-call JSON on stdin; on a policy violation it emits a PreToolUse
# "deny" decision. FAILS OPEN (allows) on any internal error so a broken hook can
# never permanently wedge committing.
#
# Scope: only fires on commits that stage model/code (code/**, audit_artifacts/**,
# or *.py). Doc-only / config-only commits pass untouched.
#
# To change the required-doc set or lift the gate: edit this file or run /hooks.

input="$(cat 2>/dev/null)"
cmd="$(printf '%s' "$input" | jq -r '.tool_input.command // ""' 2>/dev/null)"

# Only act on real `git commit` invocations (adjacent tokens — ignores `git log`,
# quoted strings, etc.). Anything else: stay silent => default allow.
printf '%s' "$cmd" | grep -Eq '(^|[;&|]+[[:space:]]*)git[[:space:]]+commit([[:space:]]|$)' || exit 0

root="$(git rev-parse --show-toplevel 2>/dev/null)" || exit 0
cd "$root" 2>/dev/null || exit 0

# Files this commit will touch: staged, plus tracked-modified if -a/--all.
staged="$(git diff --cached --name-only 2>/dev/null)"
if printf '%s' "$cmd" | grep -Eq '([[:space:]]-[[:alpha:]]*a[[:alpha:]]*([[:space:]]|$))|(--all([[:space:]]|$))'; then
  staged="$staged
$(git diff --name-only 2>/dev/null)"
fi

# Empty index (and not -a): let git produce its own "nothing to commit".
[ -n "$(printf '%s' "$staged" | tr -d '[:space:]')" ] || exit 0

# Does this commit touch model/code? If not, no doc requirement.
printf '%s\n' "$staged" | grep -Eq '^(code/|audit_artifacts/)|\.py$' || exit 0

# Which required living docs are missing from the commit?
missing=""
for d in docs/STATE.md docs/PROCESS.md docs/HANDOFF.md; do
  printf '%s\n' "$staged" | grep -qxF "$d" || missing="$missing $d"
done

[ -z "$missing" ] && exit 0

reason="Commit blocked by project policy (.claude/hooks/require-docs-before-commit.sh): model/code changes are staged, but these living docs are NOT updated in this commit:${missing}. Update each so it reflects the change — STATE.md (current state / calibration / findings), PROCESS.md (workflow / verification), HANDOFF.md (session status / priorities) — plus CLAUDE.md if instructions changed. Stage the doc updates and re-commit. To change the required set or lift this gate, edit the hook or run /hooks."

jq -cn --arg r "$reason" \
  '{hookSpecificOutput:{hookEventName:"PreToolUse",permissionDecision:"deny",permissionDecisionReason:$r}}' 2>/dev/null || exit 0
exit 0
