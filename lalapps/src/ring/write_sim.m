function write_sim(input_file,output_file) 

% read in inspiral table
  insp=readMeta(input_file,'sim_inspiral',0,'mass1,mass2,mchirp,eta,spin1x,spin1y,spin1z,spin2x,spin2y,spin2z,geocent_end_time,geocent_end_time_ns,h_end_time,h_end_time_ns,l_end_time,l_end_time_ns,end_time_gmst,distance,longitude,latitude,eff_dist_h,eff_dist_l,polarization,inclination');

% calculate the expected final spin, mass, quality, frequency using the inspiral parameters
  ring.a_final = sqrt(12).*insp.eta - 2.9.*insp.eta.^2;
  ring.M_final = ( 1 + (sqrt(8/9)-1).*insp.eta - 0.498.*insp.eta.^2 ) .* (insp.mass1+insp.mass2);
  ring.Q_ring = Qofa(ring.a_final);
  ring.f_ring = fofMa(ring.M_final,ring.a_final);

% set the following parameters to zero 
  for i=1:length(insp.mass1)
    ring.phase(i) = 0;
    ring.epsilon(i) = 0;
    ring.amp(i) = 0;
    ring.hrss(i) = 0;
    ring.hrss_h(i) = 0;
    ring.hrss_l(i) = 0;
  end

% write out the ringdown parameters to a xml file

  fid=fopen(output_file,'w');

  fprintf(fid, '<?xml version="1.0" encoding="utf-8" ?> \n');
  fprintf(fid, '<!DOCTYPE LIGO_LW SYSTEM "http://ldas-sw.ligo.caltech.edu/doc/ligolwAPI/html/ligolw_dtd.txt"><LIGO_LW>\n');

% write the sim_ringdown table: note that the times written out to the sim_ringdown table
  % are those from the sim_inspiral table.
  fprintf(fid, '<Table Name="sim_ringdowngroup:sim_ringdown:table">\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:process_id" Type="ilwd:char"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:waveform" Type="lstring"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:coordinates" Type="lstring"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:geocent_start_time" Type="int_4s"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:geocent_start_time_ns" Type="int_4s"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:h_start_time" Type="int_4s"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:h_start_time_ns" Type="int_4s"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:l_start_time" Type="int_4s"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:l_start_time_ns" Type="int_4s"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:start_time_gmst" Type="real_8"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:longitude" Type="real_4"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:latitude" Type="real_4"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:distance" Type="real_4"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:inclination" Type="real_4"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:polarization" Type="real_4"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:frequency" Type="real_4"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:quality" Type="real_4"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:phase" Type="real_4"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:mass" Type="real_4"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:spin" Type="real_4"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:epsilon" Type="real_4"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:amplitude" Type="real_4"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:eff_dist_h" Type="real_4"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:eff_dist_l" Type="real_4"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:hrss" Type="real_4"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:hrss_h" Type="real_4"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:hrss_l" Type="real_4"/>\n');
  fprintf(fid, '<Column Name="sim_ringdowngroup:sim_ringdown:simulation_id" Type="ilwd:char"/>\n');
  fprintf(fid, '<Stream Name="sim_ringdowngroup:sim_ringdown:table" Type="Local" Delimiter=",">\n');

  for i=1:length(insp.mass1)-1
    fprintf(fid, ' \t "process:process_id:0","Ringdown","EQUATORIAL",%i,%i,%i,%i,%i,%i,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,"sim_ringdown:simulation_id:%i", \n',insp.geocent_end_time(i),insp.geocent_end_time_ns(i),insp.h_end_time(i),insp.h_end_time_ns(i),insp.l_end_time(i),insp.l_end_time_ns(i),insp.end_time_gmst(i),insp.longitude(i),insp.latitude(i),insp.distance(i),insp.inclination(i),insp.polarization(i),ring.f_ring(i),ring.Q_ring(i),ring.phase(i),ring.M_final(i),ring.a_final(i),ring.epsilon(i),ring.amp(i),insp.eff_dist_h(i),insp.eff_dist_l(i),ring.hrss(i),ring.hrss_h(i),ring.hrss_l(i),i-1);
  end

% the last entry has no comma after the quotation mark
  i=i+1;
  fprintf(fid, ' \t "process:process_id:0","Ringdown","EQUATORIAL",%i,%i,%i,%i,%i,%i,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,"sim_ringdown:simulation_id:%i" \n',insp.geocent_end_time(i),insp.geocent_end_time_ns(i),insp.h_end_time(i),insp.h_end_time_ns(i),insp.l_end_time(i),insp.l_end_time_ns(i),insp.end_time_gmst(i),insp.longitude(i),insp.latitude(i),insp.distance(i),insp.inclination(i),insp.polarization(i),ring.f_ring(i),ring.Q_ring(i),ring.phase(i),ring.M_final(i),ring.a_final(i),ring.epsilon(i),ring.amp(i),insp.eff_dist_h(i),insp.eff_dist_l(i),ring.hrss(i),ring.hrss_h(i),ring.hrss_l(i),i-1);

  fprintf(fid, ' </Stream> \n');
  fprintf(fid, ' </Table> \n');
  fprintf(fid, ' </LIGO_LW> \n');
         
  fclose(fid);
%end

