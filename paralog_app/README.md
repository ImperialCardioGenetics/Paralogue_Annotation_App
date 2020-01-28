# Paralog Annotator

This repository contains code for a Shiny App to use Paralogue Variant Annotation.

A manuscript describing our approach is here (TODO:preprint-link).
You can see the app in action here (TODO:app-link).

**TODO:**

* Panda
	- Write up an ABOUT tab for the page that will link to manuscript
	- Output known ClinVar ids for query vars :heavy_check_mark:
	- Genarate url for ClinVar result variants and insert into output df as links :heavy_check_mark:
	- Generate url for ENSEMBL paralog alignemts and insert into output df as links :heavy_check_mark:
	- Fix links in download file :heavy_check_mark:
	- Fix column lenght names :heavy_check_mark:
	- Fix logo width
	- Look at changing layout to move sidepanel to top and results table below it
	- Fix Upload and Paste variant input :heavy_check_mark:
	- Error catch incorrect format in Upload and Paste inputs

* Nick
	- Swap place holder for all missense :heavy_check_mark:
	- If query variant already reported in ClinVar then return report
	- If different mutation but in same position as query variant return that
	- Code function to search for known variants in same codon
	- Return paralogous locations only - separate table tab?
	- Return Clingen gene reports - link :heavy_check_mark:
	- Query variant - return VEP CSQ? Will have to utilise paraloc output file

	
