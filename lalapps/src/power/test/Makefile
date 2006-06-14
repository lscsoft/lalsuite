all : rate_fake.png

seglist.segwizard :
	python random_segments.py >seglist.segwizard

fake.xml : seglist.segwizard
	../lalapps_power_fake --verbose --segwizard-segments seglist.segwizard >fake.xml

rate.xml : fake.xml
	../lalapps_power_trigfreq --verbose --only-rate --program lalapps_power_fake fake.xml >rate.xml

rate_fake.png : rate.xml
	../plottrigfreq --verbose rate.xml

clean :
	rm -f *.xml *.png
