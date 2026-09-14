#!/usr/bin/env bash
# UserPromptSubmit hook — standing writing-style rule for this repo.
#
# Injects a reminder that model symbols must be glossed with a short definition
# the first time they appear in prose. The model's variable names (n_inter_F,
# rdep_F, def_rate_ss, goods_mkt_D, ...) are opaque without one, and prose that
# leans on them is unreadable to anyone but the author.
#
# Emits JSON on stdout; `additionalContext` is injected into the model's context
# for the turn. Contains no double quotes so the quoted heredoc below is valid
# JSON verbatim — keep it that way when editing the text.
cat <<'JSON'
{"hookSpecificOutput":{"hookEventName":"UserPromptSubmit","additionalContext":"Writing-style rule for this repo (prose only): the first time a model variable, parameter, coefficient or residual appears in a sentence, gloss it with a short definition in parentheses immediately after it — e.g. `n_inter_F` (F-bank net worth), `rdep_F` (F ex-ante real deposit rate), `def_rate_ss` (steady-state default probability), `goods_mkt_D` (D goods-market residual), `psi_lambda_B` (collateral-friction amplification dial), `size_F` (F size relative to D). Applies to prose and table cells; skip inside code blocks, file paths, diffs and pasted program output. Do not re-gloss the same symbol twice in one reply."}}
JSON
