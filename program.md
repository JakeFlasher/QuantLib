# autoresearch

This is an experiment to have the LLM do its own research by coupling Claude Code and Codex using humanize skills.

## Setup

To perform the results analysis of our implemented approaches in QuantLib, work with the user to:

1. **Read papers in `ql/paper_references`**: Read the in-scope pdf papers and their corresponding in-scope .tex to understand what we have achieved so far.  if you find anything from the pdf or math formulations in latex ambiguous, you must use /ask-codex to enquire gpt-5.4:xhigh for detailed explanation, keep going until you're got full grasps of these papers and have no questions.
2. **Look at existing implementation**: The implementation of the proposed approaches from papers is almost done, containing a list of source files changed in `/ql` (marked with a starting tag "// r6" or similar "r6" tag) as well as in `/test-suite` (maybe marked or not marked), you can check the git commit history for a detailed list of changed files. You need to understand their logic implementation for coherence, correctness and mathematical rigor. If you found any bugs, incoherence or logical fault in the whole process, first discuss using /ask-codex, if codex also agrees with you, you MUST stop immediately and report to the user to confirm whether to fix or not.
3. **Present a plan for doing professional results analysis**: The existing test-suites as well as scripts for showing results, performance may be not sufficient for presenting in a top-tier journal. You need to use /humanize:gen-plan to generate a detailed plan that does: 1. generating vivid, professional figures for comparison; 2. writing a comprehensin markdown file consisting of comprehensive results analysis, referencing figures that are newly generated for submission to top-tier journals.
5. **Loop using /humanize:start-rlcr-loop to complete the above plan**: Evoke the skill /humanize:start-rlcr-loop and work with Codex to complete the above plan.


**What you CAN do:**
- Directly modify any files of the QuantLib source implementation directory (`/ql/*` and `/test-suite/*`) that are already modified or touched by prior implementations.
- Use all kinds of built-in skills in the /humanize: and /ask-codex. The default model for enquiring Codex is gpt-5.4:xhigh. 

**What you CANNOT do:**
- Directly modify any files of the QuantLib source implementation directory (`/ql/*` and `/test-suite/*`) that are not already modified or touched by prior implementations. If you require further interference with the source code of other parts of the QuantLib, you must ask for user approval.

**Simplicity criterion**: All else being equal, a simpler idea is better.

