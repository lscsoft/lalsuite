import sys
import glob
from optparse import OptionParser

from glue import cbcwebpage

from glue.ligolw import lsctables, table, utils
from glue.ligolw.utils import process

optp = OptionParser()
opts, args = optp.parse_args()

xmlfiles = sorted(glob.glob("*HVETO_ROUND_*.xml.gz"))

ignore_list = process.get_process_params(utils.load_filename(xmlfiles[0]), "laldetchar-hveto", "--ignore-list")[0]
reference_channel = process.get_process_params(utils.load_filename(xmlfiles[0]), "laldetchar-hveto", "--reference-channel")[0]
gps_start = process.get_process_params(utils.load_filename(xmlfiles[0]), "laldetchar-hveto", "--gps-start")[0]
gps_end = process.get_process_params(utils.load_filename(xmlfiles[0]), "laldetchar-hveto", "--gps-end")[0]

page = cbcwebpage.cbcpage(title="HVeto Summary")

subp = page.add_subpage("overview", "HVeto Run Overview", "HVeto Run Overview")

#subp.add_text(txt="GPS segment: %10.2f -- %10.2f" % (gps_start, gps_end))
subp.div("""
<big><b>
GPS segment: %10.2f -- %10.2f <br/>
Reference channel: %s <br/>
</b></big>
""" % (gps_start, gps_end, reference_channel))

imgtab = cbcwebpage.image_glob("overall_snr_hist.png")
cap = "Overall SNR histograms"
cap1 = "SNR distribution of the reference channel before and after the HVeto analysis."
subp.add_table(imgtab, cap, cap1)

imgtab = cbcwebpage.image_glob("effdead.png")
cap = "Efficiency / Deadtime"
cap1 = "Efficiency versus Deadtime after each round."
subp.add_table(imgtab, cap, cap1)

for i, xmlsum in enumerate(xmlfiles):
	i+=1

	for sv in table.get_table(utils.load_filename(xmlsum), lsctables.SearchSummVarsTable.tableName):
		if sv.name != "winner":
			continue
		winner = sv.string

	page.add_subpage("round%d" % i, "Round %d" % i, "Round %d" % i)
	sec = page.subpages["round%d" % i].add_section("r1summ", "Round %d Summary" % i)

	ttable, name = cbcwebpage.wiki_table_parse(xmlsum.replace("xml.gz", "txt"))
	sec.add_table(ttable[0], "Round summary", "Round summary")
	sec.add_link(text="Full Result XML", href=xmlsum)

	ssec = sec.add_section("r1summplot", "Round %d Summary Plot" % i)

	imgtab = cbcwebpage.image_glob("*_round_%d_summary.png" % i)
	cap = "Round %d HVeto summary" % i
	cap1 = "Round %d HVeto summary, winner %s" % (i, winner)
	ssec.add_table(imgtab, cap, cap1)

	# TODO: Add significance drop

	ssec = sec.add_section("vetotrigs", "Round %d Vetoed Triggers" % i)
	ttable, name = cbcwebpage.wiki_table_parse(xmlsum.replace("xml.gz", "trigs"))
	ssec.add_table(ttable[0], "Vetoed Triggers", "Vetoed Triggers")

	ssec = sec.add_section("vetosegs", "Round %d Veto Segments" % i)
	ttable, name = cbcwebpage.wiki_table_parse(xmlsum.replace("xml.gz", "segs"))
	ssec.add_table(ttable[0], "Veto Segments", "Veto Segments")

page.write("index")
