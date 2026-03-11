# autoresearch

This is an experiment to have the LLM do its own research by coupling Claude Code and Codex using humanize skills.

## Setup

To set up a new PhD research proposal, work with the user to:

1. **Read papers in `ql/paper_references`**: Read the in-scope pdf papers and their corresponding in-scope .tex to understand what we're trying to achieve.  if you find anything from the pdf or math formulations in latex ambiguous, you must use /ask-codex to enquire gpt-5.4:xhigh for detailed explanation, keep going until you're fully understood of these papers and have no questions. 
2. **Look at existing implementation**: The git repository is half-done already, containing a list of source files changed in `/ql` (marked with a starting tag "// r6" or similar "r6" tag). You need to first check their logic for coherence, correctness and mathematical rigor. 
3. **Discuss with Codex before making changes**: Once you have your own verdict, and adjustments to be made (if you find the current implementation wrong). Dicuss with Codex and let Codex review them first. If your approach differs from Codex or if Codex disagrees on your verdict, you must stop to ask for user to confirm who is correct.
4. **Improving the original test-suite and add more analytical/practical results**: The original implementation only touches two test files in `/test-suite`, however, the goal is to draw some nice curves and give analytical results just as in the original papers. Thus we may need to define custom scenario or parameters to showcase how the CN variant compared to the original CN in QuantLib.
5. **Loop until Codex is satisfied with your implementation**: Every time you commit any changes to the code, let Codex review it first.  If your approach differs from Codex or if Codex disagrees on your changes, you must stop to ask for user to confirm who is correct.
6. **Check papers again for extension and further future work**: List a few of future work directions such as implementing it in multi-dimensions, and discuss with Codex to review them before presenting the final roadmap to the user.
 

**What you CAN do:**
- Directly modify any files outside of the QuantLib source implementation directory (`/ql/*`). 
- Directly modify any files inside of the QuantLib source implementation directory (`/ql/*`) that has a tarting tag "// r6" or similar "r6" tag.
- Use all kinds of built-in skills in the /humanize: and /ask-codex. The default model for enquiring Codex is gpt-5.4:xhigh. \

**What you CANNOT do:**
- Directly modify any other source file implementations (`/ql/*`) of QuantLib without the starting tag "// r6" or similar "r6" tag. If you require further interference with the source code of other parts of the QuantLib, you must ask for human approval.  

**Simplicity criterion**: All else being equal, a simpler idea is better. 


