Last login: Fri Jul 11 19:24:38 on tty??
Noa:~ taniaregimbau$ ssh -l regimbau projekct.oca.eu
regimbau@projekct.oca.eu's password: 
Last login: Fri Jul 11 19:17:24 2008 from vpn-n-nat1.obs-nice.fr
Linux for High Performance Computing

================================================================
RAPPEL : L'usage des moyens informatiques hébergés à l'OCA est 
réglementé par la charte pour l'usage de ressources informatiques 
et de services internet du CNRS.
Cette charte est disponible à l'adresse :
https://sitweb.oca.eu/IMG/pdf/Charte_informatique.pdf

Pour plus de renseignements, vous pouvez consulter :
https://sitweb.oca.eu/article78.html
================================================================

[regimbau@projekct ~]$ mkdir data
[regimbau@projekct ~]$ sftp tania@norma.oca.euConnecting to norma.oca.eu...
tania@norma.oca.eu's password: 
sftp> ls
News            PUBLICATIONS    Popcorn         Popcorn.tar     Popcorn.zip     
Strings         Strings.tar     analysis        bin             codes           
codes.tar       data            documents       ldg-4.4         ligotools.tar   
lscgrid         opt             save            src             stages          
texmf           tmp             trusted.caches  vdt             
sftp> cd data
sftp> ls
PROJECT1a  p1b        
sftp> cd p1b
sftp> ls
DNS    MT     noise  
sftp> cd ..
sftp> cd ..
sftp> cd Popcorn
sftp> ls
H1.cache        H2.cache        cache.c         compil          data1           
data1.tar       data2           data2.tar       log             matlab          
output          popcorn         popcorn.c       popcorn.c_save  popcorn.c~      
sftp> cd data2
sftp> get *
Fetching /home/tania/Popcorn/data2/H2-SIM-700000000-60.gwf to H2-SIM-700000000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700000000-60 100% 3717KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700000060-60.gwf to H2-SIM-700000060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700000060-60 100% 3714KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700000120-60.gwf to H2-SIM-700000120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700000120-60 100% 3704KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700000180-60.gwf to H2-SIM-700000180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700000180-60 100% 3697KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700000240-60.gwf to H2-SIM-700000240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700000240-60 100% 3712KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700000300-60.gwf to H2-SIM-700000300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700000300-60 100% 3714KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700000360-60.gwf to H2-SIM-700000360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700000360-60 100% 3709KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700000420-60.gwf to H2-SIM-700000420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700000420-60 100% 3710KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700000480-60.gwf to H2-SIM-700000480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700000480-60 100% 3713KB 742.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700000540-60.gwf to H2-SIM-700000540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700000540-60 100% 3706KB 926.5KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700000600-60.gwf to H2-SIM-700000600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700000600-60 100% 3708KB 741.7KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700000660-60.gwf to H2-SIM-700000660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700000660-60 100% 3715KB 619.2KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700000720-60.gwf to H2-SIM-700000720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700000720-60 100% 3715KB 371.5KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700000780-60.gwf to H2-SIM-700000780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700000780-60 100% 3714KB 123.8KB/s   00:30    
Fetching /home/tania/Popcorn/data2/H2-SIM-700000840-60.gwf to H2-SIM-700000840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700000840-60 100% 3707KB 741.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700000900-60.gwf to H2-SIM-700000900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700000900-60 100% 3713KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700000960-60.gwf to H2-SIM-700000960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700000960-60 100% 3707KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700001020-60.gwf to H2-SIM-700001020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700001020-60 100% 3716KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700001080-60.gwf to H2-SIM-700001080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700001080-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700001140-60.gwf to H2-SIM-700001140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700001140-60 100% 3707KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700001200-60.gwf to H2-SIM-700001200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700001200-60 100% 3706KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700001260-60.gwf to H2-SIM-700001260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700001260-60 100% 3711KB 742.2KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700001320-60.gwf to H2-SIM-700001320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700001320-60 100% 3707KB 926.7KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700001380-60.gwf to H2-SIM-700001380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700001380-60 100% 3711KB 618.4KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700001440-60.gwf to H2-SIM-700001440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700001440-60 100% 3710KB 741.9KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700001500-60.gwf to H2-SIM-700001500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700001500-60 100% 3700KB 284.6KB/s   00:13    
Fetching /home/tania/Popcorn/data2/H2-SIM-700001560-60.gwf to H2-SIM-700001560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700001560-60 100% 3707KB 142.6KB/s   00:26    
Fetching /home/tania/Popcorn/data2/H2-SIM-700001620-60.gwf to H2-SIM-700001620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700001620-60 100% 3697KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700001680-60.gwf to H2-SIM-700001680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700001680-60 100% 3716KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700001740-60.gwf to H2-SIM-700001740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700001740-60 100% 3720KB 930.0KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700001800-60.gwf to H2-SIM-700001800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700001800-60 100% 3710KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700001860-60.gwf to H2-SIM-700001860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700001860-60 100% 3722KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700001920-60.gwf to H2-SIM-700001920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700001920-60 100% 3722KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700001980-60.gwf to H2-SIM-700001980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700001980-60 100% 3716KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700002040-60.gwf to H2-SIM-700002040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700002040-60 100% 3707KB 741.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700002100-60.gwf to H2-SIM-700002100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700002100-60 100% 3713KB 928.2KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700002160-60.gwf to H2-SIM-700002160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700002160-60 100% 3700KB 217.6KB/s   00:17    
Fetching /home/tania/Popcorn/data2/H2-SIM-700002220-60.gwf to H2-SIM-700002220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700002220-60   0%    0     0.0KB/s   --:-- ETA

/home/tania/Popcorn/data2/H2-SIM-700002220-60 100% 3716KB 177.0KB/s   00:21    
Fetching /home/tania/Popcorn/data2/H2-SIM-700002280-60.gwf to H2-SIM-700002280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700002280-60 100% 3716KB 309.6KB/s   00:12    
Fetching /home/tania/Popcorn/data2/H2-SIM-700002340-60.gwf to H2-SIM-700002340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700002340-60 100% 3707KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700002400-60.gwf to H2-SIM-700002400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700002400-60 100% 3722KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700002460-60.gwf to H2-SIM-700002460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700002460-60 100% 3718KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700002520-60.gwf to H2-SIM-700002520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700002520-60 100% 3719KB 265.6KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700002580-60.gwf to H2-SIM-700002580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700002580-60 100% 3701KB 217.7KB/s   00:17    
Fetching /home/tania/Popcorn/data2/H2-SIM-700002640-60.gwf to H2-SIM-700002640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700002640-60 100% 3727KB 745.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700002700-60.gwf to H2-SIM-700002700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700002700-60 100% 3724KB 930.9KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700002760-60.gwf to H2-SIM-700002760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700002760-60 100% 3703KB 740.6KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700002820-60.gwf to H2-SIM-700002820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700002820-60 100% 3712KB 285.5KB/s   00:13    
Fetching /home/tania/Popcorn/data2/H2-SIM-700002880-60.gwf to H2-SIM-700002880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700002880-60 100% 3708KB 185.4KB/s   00:20    
Fetching /home/tania/Popcorn/data2/H2-SIM-700002940-60.gwf to H2-SIM-700002940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700002940-60 100% 3710KB 371.0KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700003000-60.gwf to H2-SIM-700003000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700003000-60 100% 3714KB 619.0KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700003060-60.gwf to H2-SIM-700003060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700003060-60 100% 3708KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700003120-60.gwf to H2-SIM-700003120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700003120-60 100% 3706KB 926.6KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700003180-60.gwf to H2-SIM-700003180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700003180-60 100% 3714KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700003240-60.gwf to H2-SIM-700003240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700003240-60 100% 3714KB 285.7KB/s   00:13    
Fetching /home/tania/Popcorn/data2/H2-SIM-700003300-60.gwf to H2-SIM-700003300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700003300-60 100% 3709KB 195.2KB/s   00:19    
Fetching /home/tania/Popcorn/data2/H2-SIM-700003360-60.gwf to H2-SIM-700003360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700003360-60 100% 3700KB 616.6KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700003420-60.gwf to H2-SIM-700003420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700003420-60 100% 3709KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700003480-60.gwf to H2-SIM-700003480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700003480-60 100% 3714KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700003540-60.gwf to H2-SIM-700003540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700003540-60 100% 3716KB 177.0KB/s   00:21    
Fetching /home/tania/Popcorn/data2/H2-SIM-700003600-60.gwf to H2-SIM-700003600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700003600-60 100% 3706KB 195.1KB/s   00:19    
Fetching /home/tania/Popcorn/data2/H2-SIM-700003660-60.gwf to H2-SIM-700003660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700003660-60 100% 3704KB 529.1KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700003720-60.gwf to H2-SIM-700003720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700003720-60 100% 3716KB 743.1KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700003780-60.gwf to H2-SIM-700003780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700003780-60 100% 3716KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700003840-60.gwf to H2-SIM-700003840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700003840-60 100% 3715KB 265.4KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700003900-60.gwf to H2-SIM-700003900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700003900-60 100% 3701KB 185.0KB/s   00:20    
Fetching /home/tania/Popcorn/data2/H2-SIM-700003960-60.gwf to H2-SIM-700003960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700003960-60 100% 3710KB 927.5KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700004020-60.gwf to H2-SIM-700004020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700004020-60 100% 3709KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700004080-60.gwf to H2-SIM-700004080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700004080-60 100% 3704KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700004140-60.gwf to H2-SIM-700004140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700004140-60 100% 3705KB 411.6KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700004200-60.gwf to H2-SIM-700004200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700004200-60 100% 3717KB 148.7KB/s   00:25    
Fetching /home/tania/Popcorn/data2/H2-SIM-700004260-60.gwf to H2-SIM-700004260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700004260-60 100% 3689KB 527.0KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700004320-60.gwf to H2-SIM-700004320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700004320-60 100% 3719KB 743.9KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700004380-60.gwf to H2-SIM-700004380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700004380-60 100% 3707KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700004440-60.gwf to H2-SIM-700004440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700004440-60 100% 3722KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700004500-60.gwf to H2-SIM-700004500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700004500-60 100% 3699KB 246.6KB/s   00:15    
Fetching /home/tania/Popcorn/data2/H2-SIM-700004560-60.gwf to H2-SIM-700004560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700004560-60 100% 3708KB 168.6KB/s   00:22    
Fetching /home/tania/Popcorn/data2/H2-SIM-700004620-60.gwf to H2-SIM-700004620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700004620-60 100% 3714KB 928.4KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700004680-60.gwf to H2-SIM-700004680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700004680-60 100% 3718KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700004740-60.gwf to H2-SIM-700004740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700004740-60 100% 3707KB 741.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700004800-60.gwf to H2-SIM-700004800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700004800-60 100% 3711KB 265.0KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700004860-60.gwf to H2-SIM-700004860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700004860-60 100% 3700KB 119.4KB/s   00:31    
Fetching /home/tania/Popcorn/data2/H2-SIM-700004920-60.gwf to H2-SIM-700004920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700004920-60 100% 3705KB 741.1KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700004980-60.gwf to H2-SIM-700004980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700004980-60 100% 3710KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700005040-60.gwf to H2-SIM-700005040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700005040-60 100% 3715KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700005100-60.gwf to H2-SIM-700005100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700005100-60 100% 3698KB 264.1KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700005160-60.gwf to H2-SIM-700005160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700005160-60 100% 3714KB 185.7KB/s   00:20    
Fetching /home/tania/Popcorn/data2/H2-SIM-700005220-60.gwf to H2-SIM-700005220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700005220-60 100% 3700KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700005280-60.gwf to H2-SIM-700005280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700005280-60 100% 3713KB 928.3KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700005340-60.gwf to H2-SIM-700005340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700005340-60 100% 3704KB 740.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700005400-60.gwf to H2-SIM-700005400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700005400-60 100% 3707KB 926.9KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700005460-60.gwf to H2-SIM-700005460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700005460-60 100% 3704KB 264.5KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700005520-60.gwf to H2-SIM-700005520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700005520-60 100% 3703KB 176.3KB/s   00:21    
Fetching /home/tania/Popcorn/data2/H2-SIM-700005580-60.gwf to H2-SIM-700005580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700005580-60 100% 3707KB 529.5KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700005640-60.gwf to H2-SIM-700005640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700005640-60 100% 3713KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700005700-60.gwf to H2-SIM-700005700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700005700-60 100% 3718KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700005760-60.gwf to H2-SIM-700005760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700005760-60 100% 3701KB 217.7KB/s   00:17    
Fetching /home/tania/Popcorn/data2/H2-SIM-700005820-60.gwf to H2-SIM-700005820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700005820-60 100% 3698KB 616.3KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700005880-60.gwf to H2-SIM-700005880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700005880-60 100% 3708KB 741.6KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700005940-60.gwf to H2-SIM-700005940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700005940-60 100% 3715KB 412.8KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700006000-60.gwf to H2-SIM-700006000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700006000-60 100% 3718KB 154.9KB/s   00:24    
Fetching /home/tania/Popcorn/data2/H2-SIM-700006060-60.gwf to H2-SIM-700006060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700006060-60 100% 3705KB 926.3KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700006120-60.gwf to H2-SIM-700006120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700006120-60 100% 3720KB 929.9KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700006180-60.gwf to H2-SIM-700006180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700006180-60 100% 3706KB 463.3KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700006240-60.gwf to H2-SIM-700006240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700006240-60 100% 3711KB 168.7KB/s   00:22    
Fetching /home/tania/Popcorn/data2/H2-SIM-700006300-60.gwf to H2-SIM-700006300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700006300-60 100% 3716KB 743.2KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700006360-60.gwf to H2-SIM-700006360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700006360-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700006420-60.gwf to H2-SIM-700006420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700006420-60 100% 3708KB 195.2KB/s   00:19    
Fetching /home/tania/Popcorn/data2/H2-SIM-700006480-60.gwf to H2-SIM-700006480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700006480-60 100% 3703KB 264.5KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700006540-60.gwf to H2-SIM-700006540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700006540-60 100% 3716KB 371.6KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700006600-60.gwf to H2-SIM-700006600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700006600-60 100% 3704KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700006660-60.gwf to H2-SIM-700006660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700006660-60 100% 3713KB 464.1KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700006720-60.gwf to H2-SIM-700006720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700006720-60 100% 3702KB 154.3KB/s   00:24    
Fetching /home/tania/Popcorn/data2/H2-SIM-700006780-60.gwf to H2-SIM-700006780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700006780-60 100% 3703KB 740.6KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700006840-60.gwf to H2-SIM-700006840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700006840-60 100% 3708KB 741.7KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700006900-60.gwf to H2-SIM-700006900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700006900-60 100% 3714KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700006960-60.gwf to H2-SIM-700006960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700006960-60 100% 3707KB 264.8KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700007020-60.gwf to H2-SIM-700007020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700007020-60 100% 3714KB 168.8KB/s   00:22    
Fetching /home/tania/Popcorn/data2/H2-SIM-700007080-60.gwf to H2-SIM-700007080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700007080-60 100% 3720KB 531.5KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700007140-60.gwf to H2-SIM-700007140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700007140-60 100% 3706KB 741.2KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700007200-60.gwf to H2-SIM-700007200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700007200-60 100% 3711KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700007260-60.gwf to H2-SIM-700007260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700007260-60 100% 3713KB 206.3KB/s   00:18    
Fetching /home/tania/Popcorn/data2/H2-SIM-700007320-60.gwf to H2-SIM-700007320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700007320-60 100% 3711KB 185.5KB/s   00:20    
Fetching /home/tania/Popcorn/data2/H2-SIM-700007380-60.gwf to H2-SIM-700007380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700007380-60 100% 3716KB 412.9KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700007440-60.gwf to H2-SIM-700007440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700007440-60 100% 3715KB 742.9KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700007500-60.gwf to H2-SIM-700007500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700007500-60 100% 3718KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700007560-60.gwf to H2-SIM-700007560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700007560-60 100% 3698KB 739.6KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700007620-60.gwf to H2-SIM-700007620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700007620-60 100% 3717KB 265.5KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700007680-60.gwf to H2-SIM-700007680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700007680-60 100% 3701KB 284.7KB/s   00:13    
Fetching /home/tania/Popcorn/data2/H2-SIM-700007740-60.gwf to H2-SIM-700007740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700007740-60 100% 3707KB 529.5KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700007800-60.gwf to H2-SIM-700007800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700007800-60 100% 3720KB 338.2KB/s   00:11    
Fetching /home/tania/Popcorn/data2/H2-SIM-700007860-60.gwf to H2-SIM-700007860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700007860-60 100% 3699KB 194.7KB/s   00:19    
Fetching /home/tania/Popcorn/data2/H2-SIM-700007920-60.gwf to H2-SIM-700007920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700007920-60 100% 3716KB 371.6KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700007980-60.gwf to H2-SIM-700007980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700007980-60 100% 3716KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700008040-60.gwf to H2-SIM-700008040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700008040-60 100% 3704KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700008100-60.gwf to H2-SIM-700008100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700008100-60 100% 3712KB 265.1KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700008160-60.gwf to H2-SIM-700008160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700008160-60 100% 3709KB 206.0KB/s   00:18    
Fetching /home/tania/Popcorn/data2/H2-SIM-700008220-60.gwf to H2-SIM-700008220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700008220-60 100% 3708KB 529.8KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700008280-60.gwf to H2-SIM-700008280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700008280-60 100% 3690KB 738.0KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700008340-60.gwf to H2-SIM-700008340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700008340-60 100% 3714KB 742.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700008400-60.gwf to H2-SIM-700008400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700008400-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700008460-60.gwf to H2-SIM-700008460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700008460-60 100% 3716KB 265.4KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700008520-60.gwf to H2-SIM-700008520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700008520-60 100% 3724KB 137.9KB/s   00:27    
Fetching /home/tania/Popcorn/data2/H2-SIM-700008580-60.gwf to H2-SIM-700008580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700008580-60 100% 3707KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700008640-60.gwf to H2-SIM-700008640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700008640-60 100% 3713KB 742.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700008700-60.gwf to H2-SIM-700008700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700008700-60 100% 3693KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700008760-60.gwf to H2-SIM-700008760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700008760-60 100% 3697KB 264.1KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700008820-60.gwf to H2-SIM-700008820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700008820-60 100% 3722KB 177.3KB/s   00:21    
Fetching /home/tania/Popcorn/data2/H2-SIM-700008880-60.gwf to H2-SIM-700008880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700008880-60 100% 3701KB 616.8KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700008940-60.gwf to H2-SIM-700008940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700008940-60 100% 3708KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700009000-60.gwf to H2-SIM-700009000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700009000-60 100% 3702KB 740.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700009060-60.gwf to H2-SIM-700009060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700009060-60 100% 3708KB 195.2KB/s   00:19    
Fetching /home/tania/Popcorn/data2/H2-SIM-700009120-60.gwf to H2-SIM-700009120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700009120-60 100% 3709KB 463.6KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700009180-60.gwf to H2-SIM-700009180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700009180-60 100% 3717KB 247.8KB/s   00:15    
Fetching /home/tania/Popcorn/data2/H2-SIM-700009240-60.gwf to H2-SIM-700009240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700009240-60 100% 3705KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700009300-60.gwf to H2-SIM-700009300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700009300-60 100% 3705KB 463.2KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700009360-60.gwf to H2-SIM-700009360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700009360-60 100% 3701KB 168.2KB/s   00:22    
Fetching /home/tania/Popcorn/data2/H2-SIM-700009420-60.gwf to H2-SIM-700009420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700009420-60 100% 3715KB 743.0KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700009480-60.gwf to H2-SIM-700009480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700009480-60 100% 3718KB 619.7KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700009540-60.gwf to H2-SIM-700009540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700009540-60 100% 3706KB 247.1KB/s   00:15    
Fetching /home/tania/Popcorn/data2/H2-SIM-700009600-60.gwf to H2-SIM-700009600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700009600-60 100% 3709KB 206.1KB/s   00:18    
Fetching /home/tania/Popcorn/data2/H2-SIM-700009660-60.gwf to H2-SIM-700009660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700009660-60 100% 3714KB 742.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700009720-60.gwf to H2-SIM-700009720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700009720-60 100% 3707KB 411.9KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700009780-60.gwf to H2-SIM-700009780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700009780-60 100% 3724KB 186.2KB/s   00:20    
Fetching /home/tania/Popcorn/data2/H2-SIM-700009840-60.gwf to H2-SIM-700009840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700009840-60 100% 3713KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700009900-60.gwf to H2-SIM-700009900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700009900-60 100% 3713KB 742.7KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700009960-60.gwf to H2-SIM-700009960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700009960-60 100% 3699KB 462.4KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700010020-60.gwf to H2-SIM-700010020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700010020-60   0%    0     0.0KB/s   --:-- ETA




/home/tania/Popcorn/data2/H2-SIM-700010020-60 100% 3707KB 195.1KB/s   00:19    
Fetching /home/tania/Popcorn/data2/H2-SIM-700010080-60.gwf to H2-SIM-700010080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700010080-60 100% 3707KB 926.8KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700010140-60.gwf to H2-SIM-700010140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700010140-60 100% 3701KB 740.1KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700010200-60.gwf to H2-SIM-700010200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700010200-60 100% 3705KB 264.6KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700010260-60.gwf to H2-SIM-700010260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700010260-60 100% 3703KB 161.0KB/s   00:23    
Fetching /home/tania/Popcorn/data2/H2-SIM-700010320-60.gwf to H2-SIM-700010320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700010320-60 100% 3712KB 618.7KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700010380-60.gwf to H2-SIM-700010380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700010380-60 100% 3716KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700010440-60.gwf to H2-SIM-700010440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700010440-60 100% 3713KB 742.7KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700010500-60.gwf to H2-SIM-700010500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700010500-60 100% 3705KB 217.9KB/s   00:17    
Fetching /home/tania/Popcorn/data2/H2-SIM-700010560-60.gwf to H2-SIM-700010560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700010560-60 100% 3708KB 741.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700010620-60.gwf to H2-SIM-700010620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700010620-60 100% 3708KB 741.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700010680-60.gwf to H2-SIM-700010680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700010680-60 100% 3708KB 309.0KB/s   00:12    
Fetching /home/tania/Popcorn/data2/H2-SIM-700010740-60.gwf to H2-SIM-700010740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700010740-60 100% 3699KB 411.0KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700010800-60.gwf to H2-SIM-700010800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700010800-60 100% 3707KB 370.7KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700010860-60.gwf to H2-SIM-700010860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700010860-60 100% 3720KB 155.0KB/s   00:24    
Fetching /home/tania/Popcorn/data2/H2-SIM-700010920-60.gwf to H2-SIM-700010920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700010920-60 100% 3704KB 740.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700010980-60.gwf to H2-SIM-700010980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700010980-60 100% 3708KB 617.9KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700011040-60.gwf to H2-SIM-700011040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700011040-60 100% 3718KB 464.7KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700011100-60.gwf to H2-SIM-700011100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700011100-60 100% 3701KB 142.4KB/s   00:26    
Fetching /home/tania/Popcorn/data2/H2-SIM-700011160-60.gwf to H2-SIM-700011160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700011160-60 100% 3715KB 619.2KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700011220-60.gwf to H2-SIM-700011220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700011220-60 100% 3703KB 740.7KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700011280-60.gwf to H2-SIM-700011280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700011280-60 100% 3718KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700011340-60.gwf to H2-SIM-700011340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700011340-60 100% 3712KB 265.1KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700011400-60.gwf to H2-SIM-700011400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700011400-60 100% 3711KB 185.6KB/s   00:20    
Fetching /home/tania/Popcorn/data2/H2-SIM-700011460-60.gwf to H2-SIM-700011460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700011460-60 100% 3712KB 412.4KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700011520-60.gwf to H2-SIM-700011520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700011520-60 100% 3698KB 739.7KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700011580-60.gwf to H2-SIM-700011580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700011580-60 100% 3706KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700011640-60.gwf to H2-SIM-700011640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700011640-60 100% 3720KB 248.0KB/s   00:15    
Fetching /home/tania/Popcorn/data2/H2-SIM-700011700-60.gwf to H2-SIM-700011700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700011700-60 100% 3712KB 176.8KB/s   00:21    
Fetching /home/tania/Popcorn/data2/H2-SIM-700011760-60.gwf to H2-SIM-700011760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700011760-60 100% 3709KB 412.1KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700011820-60.gwf to H2-SIM-700011820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700011820-60 100% 3717KB 929.3KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700011880-60.gwf to H2-SIM-700011880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700011880-60 100% 3709KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700011940-60.gwf to H2-SIM-700011940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700011940-60 100% 3713KB 265.2KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700012000-60.gwf to H2-SIM-700012000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700012000-60 100% 3721KB 195.8KB/s   00:19    
Fetching /home/tania/Popcorn/data2/H2-SIM-700012060-60.gwf to H2-SIM-700012060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700012060-60 100% 3718KB 413.2KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700012120-60.gwf to H2-SIM-700012120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700012120-60 100% 3711KB 742.2KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700012180-60.gwf to H2-SIM-700012180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700012180-60 100% 3712KB 928.1KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700012240-60.gwf to H2-SIM-700012240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700012240-60 100% 3704KB 285.0KB/s   00:13    
Fetching /home/tania/Popcorn/data2/H2-SIM-700012300-60.gwf to H2-SIM-700012300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700012300-60 100% 3703KB 176.3KB/s   00:21    
Fetching /home/tania/Popcorn/data2/H2-SIM-700012360-60.gwf to H2-SIM-700012360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700012360-60 100% 3706KB 463.3KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700012420-60.gwf to H2-SIM-700012420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700012420-60 100% 3702KB 740.4KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700012480-60.gwf to H2-SIM-700012480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700012480-60 100% 3712KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700012540-60.gwf to H2-SIM-700012540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700012540-60 100% 3714KB 218.5KB/s   00:17    
Fetching /home/tania/Popcorn/data2/H2-SIM-700012600-60.gwf to H2-SIM-700012600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700012600-60 100% 3716KB 177.0KB/s   00:21    
Fetching /home/tania/Popcorn/data2/H2-SIM-700012660-60.gwf to H2-SIM-700012660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700012660-60 100% 3711KB 463.9KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700012720-60.gwf to H2-SIM-700012720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700012720-60 100% 3708KB 618.0KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700012780-60.gwf to H2-SIM-700012780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700012780-60 100% 3716KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700012840-60.gwf to H2-SIM-700012840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700012840-60 100% 3721KB 372.1KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700012900-60.gwf to H2-SIM-700012900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700012900-60 100% 3706KB 185.3KB/s   00:20    
Fetching /home/tania/Popcorn/data2/H2-SIM-700012960-60.gwf to H2-SIM-700012960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700012960-60 100% 3706KB 529.5KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700013020-60.gwf to H2-SIM-700013020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700013020-60 100% 3698KB 616.3KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700013080-60.gwf to H2-SIM-700013080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700013080-60 100% 3703KB 740.6KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700013140-60.gwf to H2-SIM-700013140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700013140-60 100% 3704KB 246.9KB/s   00:15    
Fetching /home/tania/Popcorn/data2/H2-SIM-700013200-60.gwf to H2-SIM-700013200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700013200-60 100% 3698KB 168.1KB/s   00:22    
Fetching /home/tania/Popcorn/data2/H2-SIM-700013260-60.gwf to H2-SIM-700013260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700013260-60 100% 3704KB 463.1KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700013320-60.gwf to H2-SIM-700013320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700013320-60 100% 3714KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700013380-60.gwf to H2-SIM-700013380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700013380-60 100% 3711KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700013440-60.gwf to H2-SIM-700013440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700013440-60 100% 3703KB 336.7KB/s   00:11    
Fetching /home/tania/Popcorn/data2/H2-SIM-700013500-60.gwf to H2-SIM-700013500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700013500-60 100% 3716KB 195.6KB/s   00:19    
Fetching /home/tania/Popcorn/data2/H2-SIM-700013560-60.gwf to H2-SIM-700013560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700013560-60 100% 3705KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700013620-60.gwf to H2-SIM-700013620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700013620-60 100% 3707KB 926.7KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700013680-60.gwf to H2-SIM-700013680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700013680-60 100% 3708KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700013740-60.gwf to H2-SIM-700013740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700013740-60 100% 3701KB 740.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700013800-60.gwf to H2-SIM-700013800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700013800-60 100% 3709KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700013860-60.gwf to H2-SIM-700013860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700013860-60 100% 3710KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700013920-60.gwf to H2-SIM-700013920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700013920-60 100% 3708KB 741.6KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700013980-60.gwf to H2-SIM-700013980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700013980-60 100% 3720KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700014040-60.gwf to H2-SIM-700014040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700014040-60 100% 3723KB 372.3KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700014100-60.gwf to H2-SIM-700014100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700014100-60 100% 3708KB 115.9KB/s   00:32    
Fetching /home/tania/Popcorn/data2/H2-SIM-700014160-60.gwf to H2-SIM-700014160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700014160-60 100% 3716KB 743.2KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700014220-60.gwf to H2-SIM-700014220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700014220-60 100% 3716KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700014280-60.gwf to H2-SIM-700014280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700014280-60 100% 3709KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700014340-60.gwf to H2-SIM-700014340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700014340-60 100% 3703KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700014400-60.gwf to H2-SIM-700014400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700014400-60 100% 3699KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700014460-60.gwf to H2-SIM-700014460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700014460-60 100% 3693KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700014520-60.gwf to H2-SIM-700014520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700014520-60 100% 3719KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700014580-60.gwf to H2-SIM-700014580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700014580-60 100% 3699KB 739.7KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700014640-60.gwf to H2-SIM-700014640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700014640-60 100% 3713KB 928.3KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700014700-60.gwf to H2-SIM-700014700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700014700-60 100% 3710KB 463.7KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700014760-60.gwf to H2-SIM-700014760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700014760-60 100% 3714KB 142.9KB/s   00:26    
Fetching /home/tania/Popcorn/data2/H2-SIM-700014820-60.gwf to H2-SIM-700014820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700014820-60 100% 3716KB 743.2KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700014880-60.gwf to H2-SIM-700014880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700014880-60 100% 3712KB 618.7KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700014940-60.gwf to H2-SIM-700014940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700014940-60 100% 3713KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700015000-60.gwf to H2-SIM-700015000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700015000-60 100% 3701KB 264.4KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700015060-60.gwf to H2-SIM-700015060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700015060-60 100% 3694KB 217.3KB/s   00:17    
Fetching /home/tania/Popcorn/data2/H2-SIM-700015120-60.gwf to H2-SIM-700015120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700015120-60 100% 3706KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700015180-60.gwf to H2-SIM-700015180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700015180-60 100% 3714KB 928.5KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700015240-60.gwf to H2-SIM-700015240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700015240-60 100% 3704KB 308.7KB/s   00:12    
Fetching /home/tania/Popcorn/data2/H2-SIM-700015300-60.gwf to H2-SIM-700015300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700015300-60 100% 3715KB 195.5KB/s   00:19    
Fetching /home/tania/Popcorn/data2/H2-SIM-700015360-60.gwf to H2-SIM-700015360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700015360-60 100% 3709KB 370.9KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700015420-60.gwf to H2-SIM-700015420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700015420-60 100% 3709KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700015480-60.gwf to H2-SIM-700015480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700015480-60 100% 3712KB 161.4KB/s   00:23    
Fetching /home/tania/Popcorn/data2/H2-SIM-700015540-60.gwf to H2-SIM-700015540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700015540-60 100% 3710KB 195.3KB/s   00:19    
Fetching /home/tania/Popcorn/data2/H2-SIM-700015600-60.gwf to H2-SIM-700015600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700015600-60 100% 3700KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700015660-60.gwf to H2-SIM-700015660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700015660-60 100% 3709KB 618.2KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700015720-60.gwf to H2-SIM-700015720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700015720-60 100% 3711KB 927.7KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700015780-60.gwf to H2-SIM-700015780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700015780-60 100% 3715KB 530.7KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700015840-60.gwf to H2-SIM-700015840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700015840-60 100% 3718KB 169.0KB/s   00:22    
Fetching /home/tania/Popcorn/data2/H2-SIM-700015900-60.gwf to H2-SIM-700015900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700015900-60 100% 3711KB 742.2KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700015960-60.gwf to H2-SIM-700015960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700015960-60 100% 3711KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700016020-60.gwf to H2-SIM-700016020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700016020-60 100% 3725KB 248.4KB/s   00:15    
Fetching /home/tania/Popcorn/data2/H2-SIM-700016080-60.gwf to H2-SIM-700016080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700016080-60 100% 3710KB 371.0KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700016140-60.gwf to H2-SIM-700016140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700016140-60 100% 3707KB 337.0KB/s   00:11    
Fetching /home/tania/Popcorn/data2/H2-SIM-700016200-60.gwf to H2-SIM-700016200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700016200-60 100% 3710KB 218.2KB/s   00:17    
Fetching /home/tania/Popcorn/data2/H2-SIM-700016260-60.gwf to H2-SIM-700016260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700016260-60 100% 3709KB 741.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700016320-60.gwf to H2-SIM-700016320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700016320-60 100% 3709KB 741.9KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700016380-60.gwf to H2-SIM-700016380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700016380-60 100% 3705KB 195.0KB/s   00:19    
Fetching /home/tania/Popcorn/data2/H2-SIM-700016440-60.gwf to H2-SIM-700016440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700016440-60 100% 3711KB 176.7KB/s   00:21    
Fetching /home/tania/Popcorn/data2/H2-SIM-700016500-60.gwf to H2-SIM-700016500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700016500-60 100% 3700KB 616.7KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700016560-60.gwf to H2-SIM-700016560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700016560-60 100% 3698KB 739.7KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700016620-60.gwf to H2-SIM-700016620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700016620-60 100% 3715KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700016680-60.gwf to H2-SIM-700016680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700016680-60 100% 3708KB 195.2KB/s   00:19    
Fetching /home/tania/Popcorn/data2/H2-SIM-700016740-60.gwf to H2-SIM-700016740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700016740-60 100% 3719KB 531.3KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700016800-60.gwf to H2-SIM-700016800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700016800-60 100% 3706KB 741.2KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700016860-60.gwf to H2-SIM-700016860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700016860-60 100% 3707KB 218.1KB/s   00:17    
Fetching /home/tania/Popcorn/data2/H2-SIM-700016920-60.gwf to H2-SIM-700016920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700016920-60 100% 3698KB 246.5KB/s   00:15    
Fetching /home/tania/Popcorn/data2/H2-SIM-700016980-60.gwf to H2-SIM-700016980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700016980-60 100% 3712KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700017040-60.gwf to H2-SIM-700017040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700017040-60 100% 3721KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700017100-60.gwf to H2-SIM-700017100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700017100-60 100% 3715KB 530.7KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700017160-60.gwf to H2-SIM-700017160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700017160-60 100% 3714KB 154.8KB/s   00:24    
Fetching /home/tania/Popcorn/data2/H2-SIM-700017220-60.gwf to H2-SIM-700017220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700017220-60 100% 3701KB 925.4KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700017280-60.gwf to H2-SIM-700017280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700017280-60 100% 3724KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700017340-60.gwf to H2-SIM-700017340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700017340-60 100% 3713KB 285.7KB/s   00:13    
Fetching /home/tania/Popcorn/data2/H2-SIM-700017400-60.gwf to H2-SIM-700017400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700017400-60 100% 3704KB 137.2KB/s   00:27    
Fetching /home/tania/Popcorn/data2/H2-SIM-700017460-60.gwf to H2-SIM-700017460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700017460-60 100% 3711KB 618.5KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700017520-60.gwf to H2-SIM-700017520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700017520-60 100% 3706KB 741.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700017580-60.gwf to H2-SIM-700017580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700017580-60 100% 3695KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700017640-60.gwf to H2-SIM-700017640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700017640-60 100% 3714KB 232.1KB/s   00:16    
Fetching /home/tania/Popcorn/data2/H2-SIM-700017700-60.gwf to H2-SIM-700017700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700017700-60 100% 3702KB 176.3KB/s   00:21    
Fetching /home/tania/Popcorn/data2/H2-SIM-700017760-60.gwf to H2-SIM-700017760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700017760-60 100% 3703KB 370.3KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700017820-60.gwf to H2-SIM-700017820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700017820-60 100% 3705KB 926.4KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700017880-60.gwf to H2-SIM-700017880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700017880-60 100% 3702KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700017940-60.gwf to H2-SIM-700017940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700017940-60 100% 3705KB 217.9KB/s   00:17    
Fetching /home/tania/Popcorn/data2/H2-SIM-700018000-60.gwf to H2-SIM-700018000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700018000-60 100% 3709KB 168.6KB/s   00:22    
Fetching /home/tania/Popcorn/data2/H2-SIM-700018060-60.gwf to H2-SIM-700018060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700018060-60 100% 3698KB 528.3KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700018120-60.gwf to H2-SIM-700018120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700018120-60 100% 3715KB 743.0KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700018180-60.gwf to H2-SIM-700018180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700018180-60 100% 3696KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700018240-60.gwf to H2-SIM-700018240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700018240-60 100% 3705KB 926.1KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700018300-60.gwf to H2-SIM-700018300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700018300-60 100% 3717KB 123.9KB/s   00:30    
Fetching /home/tania/Popcorn/data2/H2-SIM-700018360-60.gwf to H2-SIM-700018360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700018360-60 100% 3719KB 413.2KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700018420-60.gwf to H2-SIM-700018420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700018420-60 100% 3697KB 924.3KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700018480-60.gwf to H2-SIM-700018480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700018480-60 100% 3706KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700018540-60.gwf to H2-SIM-700018540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700018540-60 100% 3720KB 206.7KB/s   00:18    
Fetching /home/tania/Popcorn/data2/H2-SIM-700018600-60.gwf to H2-SIM-700018600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700018600-60 100% 3704KB 185.2KB/s   00:20    
Fetching /home/tania/Popcorn/data2/H2-SIM-700018660-60.gwf to H2-SIM-700018660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700018660-60 100% 3722KB 620.4KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700018720-60.gwf to H2-SIM-700018720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700018720-60 100% 3707KB 617.9KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700018780-60.gwf to H2-SIM-700018780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700018780-60 100% 3713KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700018840-60.gwf to H2-SIM-700018840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700018840-60 100% 3694KB 246.3KB/s   00:15    
Fetching /home/tania/Popcorn/data2/H2-SIM-700018900-60.gwf to H2-SIM-700018900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700018900-60 100% 3704KB 154.3KB/s   00:24    
Fetching /home/tania/Popcorn/data2/H2-SIM-700018960-60.gwf to H2-SIM-700018960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700018960-60 100% 3723KB 465.3KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700019020-60.gwf to H2-SIM-700019020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700019020-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700019080-60.gwf to H2-SIM-700019080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700019080-60 100% 3705KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700019140-60.gwf to H2-SIM-700019140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700019140-60 100% 3713KB 265.2KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700019200-60.gwf to H2-SIM-700019200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700019200-60   0%    0     0.0KB/s   --:-- ETA


/home/tania/Popcorn/data2/H2-SIM-700019200-60 100% 3704KB 176.4KB/s   00:21    
Fetching /home/tania/Popcorn/data2/H2-SIM-700019260-60.gwf to H2-SIM-700019260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700019260-60 100% 3706KB 336.9KB/s   00:11    
Fetching /home/tania/Popcorn/data2/H2-SIM-700019320-60.gwf to H2-SIM-700019320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700019320-60 100% 3710KB 741.9KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700019380-60.gwf to H2-SIM-700019380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700019380-60 100% 3699KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700019440-60.gwf to H2-SIM-700019440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700019440-60 100% 3716KB 218.6KB/s   00:17    
Fetching /home/tania/Popcorn/data2/H2-SIM-700019500-60.gwf to H2-SIM-700019500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700019500-60 100% 3707KB 168.5KB/s   00:22    
Fetching /home/tania/Popcorn/data2/H2-SIM-700019560-60.gwf to H2-SIM-700019560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700019560-60 100% 3709KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700019620-60.gwf to H2-SIM-700019620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700019620-60 100% 3711KB 742.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700019680-60.gwf to H2-SIM-700019680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700019680-60 100% 3715KB 928.8KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700019740-60.gwf to H2-SIM-700019740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700019740-60 100% 3712KB 142.8KB/s   00:26    
Fetching /home/tania/Popcorn/data2/H2-SIM-700019800-60.gwf to H2-SIM-700019800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700019800-60 100% 3712KB 371.2KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700019860-60.gwf to H2-SIM-700019860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700019860-60 100% 3711KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700019920-60.gwf to H2-SIM-700019920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700019920-60 100% 3714KB 530.5KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700019980-60.gwf to H2-SIM-700019980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700019980-60 100% 3709KB 741.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700020040-60.gwf to H2-SIM-700020040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700020040-60 100% 3703KB 194.9KB/s   00:19    
Fetching /home/tania/Popcorn/data2/H2-SIM-700020100-60.gwf to H2-SIM-700020100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700020100-60 100% 3713KB 265.2KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700020160-60.gwf to H2-SIM-700020160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700020160-60 100% 3697KB 369.7KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700020220-60.gwf to H2-SIM-700020220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700020220-60 100% 3706KB 617.7KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700020280-60.gwf to H2-SIM-700020280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700020280-60 100% 3723KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700020340-60.gwf to H2-SIM-700020340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700020340-60 100% 3704KB 308.7KB/s   00:12    
Fetching /home/tania/Popcorn/data2/H2-SIM-700020400-60.gwf to H2-SIM-700020400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700020400-60 100% 3711KB 371.1KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700020460-60.gwf to H2-SIM-700020460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700020460-60 100% 3700KB 462.4KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700020520-60.gwf to H2-SIM-700020520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700020520-60 100% 3713KB 137.5KB/s   00:27    
Fetching /home/tania/Popcorn/data2/H2-SIM-700020580-60.gwf to H2-SIM-700020580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700020580-60 100% 3700KB 616.8KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700020640-60.gwf to H2-SIM-700020640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700020640-60 100% 3697KB 739.4KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700020700-60.gwf to H2-SIM-700020700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700020700-60 100% 3712KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700020760-60.gwf to H2-SIM-700020760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700020760-60 100% 3711KB 185.6KB/s   00:20    
Fetching /home/tania/Popcorn/data2/H2-SIM-700020820-60.gwf to H2-SIM-700020820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700020820-60 100% 3706KB 370.6KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700020880-60.gwf to H2-SIM-700020880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700020880-60 100% 3712KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700020940-60.gwf to H2-SIM-700020940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700020940-60 100% 3701KB 370.1KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700021000-60.gwf to H2-SIM-700021000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700021000-60 100% 3717KB 137.7KB/s   00:27    
Fetching /home/tania/Popcorn/data2/H2-SIM-700021060-60.gwf to H2-SIM-700021060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700021060-60 100% 3708KB 927.0KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700021120-60.gwf to H2-SIM-700021120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700021120-60 100% 3696KB 528.0KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700021180-60.gwf to H2-SIM-700021180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700021180-60 100% 3709KB 264.9KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700021240-60.gwf to H2-SIM-700021240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700021240-60 100% 3710KB 247.3KB/s   00:15    
Fetching /home/tania/Popcorn/data2/H2-SIM-700021300-60.gwf to H2-SIM-700021300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700021300-60 100% 3722KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700021360-60.gwf to H2-SIM-700021360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700021360-60 100% 3705KB 168.4KB/s   00:22    
Fetching /home/tania/Popcorn/data2/H2-SIM-700021420-60.gwf to H2-SIM-700021420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700021420-60   0%    0     0.0KB/s   --:-- ETA


/home/tania/Popcorn/data2/H2-SIM-700021420-60 100% 3704KB 142.5KB/s   00:26    
Fetching /home/tania/Popcorn/data2/H2-SIM-700021480-60.gwf to H2-SIM-700021480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700021480-60 100% 3711KB 742.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700021540-60.gwf to H2-SIM-700021540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700021540-60 100% 3720KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700021600-60.gwf to H2-SIM-700021600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700021600-60 100% 3710KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700021660-60.gwf to H2-SIM-700021660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700021660-60 100% 3704KB 247.0KB/s   00:15    
Fetching /home/tania/Popcorn/data2/H2-SIM-700021720-60.gwf to H2-SIM-700021720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700021720-60 100% 3711KB 185.6KB/s   00:20    
Fetching /home/tania/Popcorn/data2/H2-SIM-700021780-60.gwf to H2-SIM-700021780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700021780-60 100% 3702KB 462.7KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700021840-60.gwf to H2-SIM-700021840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700021840-60 100% 3701KB 740.2KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700021900-60.gwf to H2-SIM-700021900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700021900-60 100% 3703KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700021960-60.gwf to H2-SIM-700021960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700021960-60 100% 3700KB 185.0KB/s   00:20    
Fetching /home/tania/Popcorn/data2/H2-SIM-700022020-60.gwf to H2-SIM-700022020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700022020-60 100% 3720KB 232.5KB/s   00:16    
Fetching /home/tania/Popcorn/data2/H2-SIM-700022080-60.gwf to H2-SIM-700022080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700022080-60 100% 3712KB 337.4KB/s   00:11    
Fetching /home/tania/Popcorn/data2/H2-SIM-700022140-60.gwf to H2-SIM-700022140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700022140-60 100% 3711KB 927.8KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700022200-60.gwf to H2-SIM-700022200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700022200-60 100% 3712KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700022260-60.gwf to H2-SIM-700022260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700022260-60 100% 3707KB 142.6KB/s   00:26    
Fetching /home/tania/Popcorn/data2/H2-SIM-700022320-60.gwf to H2-SIM-700022320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700022320-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700022320-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700022380-60.gwf to H2-SIM-700022380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700022380-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700022380-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700022440-60.gwf to H2-SIM-700022440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700022440-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700022440-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700022500-60.gwf to H2-SIM-700022500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700022500-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700022500-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700022560-60.gwf to H2-SIM-700022560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700022560-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700022560-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700022620-60.gwf to H2-SIM-700022620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700022620-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700022620-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700022680-60.gwf to H2-SIM-700022680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700022680-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700022680-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700022740-60.gwf to H2-SIM-700022740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700022740-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700022740-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700022800-60.gwf to H2-SIM-700022800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700022800-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700022800-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700022860-60.gwf to H2-SIM-700022860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700022860-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700022860-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700022920-60.gwf to H2-SIM-700022920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700022920-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700022920-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700022980-60.gwf to H2-SIM-700022980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700022980-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700022980-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700023040-60.gwf to H2-SIM-700023040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700023040-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700023040-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700023100-60.gwf to H2-SIM-700023100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700023100-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700023100-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700023160-60.gwf to H2-SIM-700023160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700023160-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700023160-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700023220-60.gwf to H2-SIM-700023220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700023220-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700023220-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700023280-60.gwf to H2-SIM-700023280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700023280-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700023280-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700023340-60.gwf to H2-SIM-700023340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700023340-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700023340-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700023400-60.gwf to H2-SIM-700023400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700023400-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700023400-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700023460-60.gwf to H2-SIM-700023460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700023460-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700023460-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700023520-60.gwf to H2-SIM-700023520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700023520-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700023520-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700023580-60.gwf to H2-SIM-700023580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700023580-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700023580-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700023640-60.gwf to H2-SIM-700023640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700023640-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700023640-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700023700-60.gwf to H2-SIM-700023700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700023700-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700023700-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700023760-60.gwf to H2-SIM-700023760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700023760-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700023760-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700023820-60.gwf to H2-SIM-700023820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700023820-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700023820-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700023880-60.gwf to H2-SIM-700023880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700023880-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700023880-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700023940-60.gwf to H2-SIM-700023940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700023940-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700023940-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700024000-60.gwf to H2-SIM-700024000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700024000-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700024000-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700024060-60.gwf to H2-SIM-700024060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700024060-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700024060-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700024120-60.gwf to H2-SIM-700024120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700024120-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700024120-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700024180-60.gwf to H2-SIM-700024180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700024180-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700024180-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700024240-60.gwf to H2-SIM-700024240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700024240-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700024240-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700024300-60.gwf to H2-SIM-700024300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700024300-60 100% 3698KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700024360-60.gwf to H2-SIM-700024360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700024360-60 100% 3712KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700024420-60.gwf to H2-SIM-700024420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700024420-60 100% 3719KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700024480-60.gwf to H2-SIM-700024480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700024480-60 100% 3710KB 463.8KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700024540-60.gwf to H2-SIM-700024540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700024540-60 100% 3707KB 411.9KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700024600-60.gwf to H2-SIM-700024600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700024600-60 100% 3717KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700024660-60.gwf to H2-SIM-700024660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700024660-60 100% 3704KB 529.1KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700024720-60.gwf to H2-SIM-700024720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700024720-60 100% 3703KB 132.3KB/s   00:28    
Fetching /home/tania/Popcorn/data2/H2-SIM-700024780-60.gwf to H2-SIM-700024780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700024780-60 100% 3702KB 740.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700024840-60.gwf to H2-SIM-700024840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700024840-60 100% 3703KB 740.6KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700024900-60.gwf to H2-SIM-700024900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700024900-60 100% 3709KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700024960-60.gwf to H2-SIM-700024960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700024960-60 100% 3711KB 927.8KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700025020-60.gwf to H2-SIM-700025020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700025020-60 100% 3702KB 925.5KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700025080-60.gwf to H2-SIM-700025080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700025080-60 100% 3704KB 926.1KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700025140-60.gwf to H2-SIM-700025140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700025140-60 100% 3699KB 284.5KB/s   00:13    
Fetching /home/tania/Popcorn/data2/H2-SIM-700025200-60.gwf to H2-SIM-700025200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700025200-60 100% 3705KB 926.2KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700025260-60.gwf to H2-SIM-700025260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700025260-60 100% 3705KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700025320-60.gwf to H2-SIM-700025320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700025320-60 100% 3716KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700025380-60.gwf to H2-SIM-700025380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700025380-60 100% 3711KB 742.2KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700025440-60.gwf to H2-SIM-700025440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700025440-60 100% 3712KB 206.2KB/s   00:18    
Fetching /home/tania/Popcorn/data2/H2-SIM-700025500-60.gwf to H2-SIM-700025500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700025500-60 100% 3714KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700025560-60.gwf to H2-SIM-700025560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700025560-60 100% 3719KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700025620-60.gwf to H2-SIM-700025620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700025620-60 100% 3712KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700025680-60.gwf to H2-SIM-700025680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700025680-60 100% 3709KB 618.1KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700025740-60.gwf to H2-SIM-700025740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700025740-60 100% 3717KB 137.7KB/s   00:27    
Fetching /home/tania/Popcorn/data2/H2-SIM-700025800-60.gwf to H2-SIM-700025800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700025800-60 100% 3707KB 741.4KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700025860-60.gwf to H2-SIM-700025860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700025860-60 100% 3710KB 927.6KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700025920-60.gwf to H2-SIM-700025920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700025920-60 100% 3708KB 927.0KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700025980-60.gwf to H2-SIM-700025980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700025980-60 100% 3714KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700026040-60.gwf to H2-SIM-700026040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700026040-60 100% 3716KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700026100-60.gwf to H2-SIM-700026100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700026100-60 100% 3706KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700026160-60.gwf to H2-SIM-700026160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700026160-60 100% 3703KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700026220-60.gwf to H2-SIM-700026220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700026220-60 100% 3709KB 741.9KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700026280-60.gwf to H2-SIM-700026280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700026280-60 100% 3700KB 925.0KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700026340-60.gwf to H2-SIM-700026340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700026340-60 100% 3712KB 928.1KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700026400-60.gwf to H2-SIM-700026400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700026400-60 100% 3711KB 742.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700026460-60.gwf to H2-SIM-700026460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700026460-60 100% 3717KB 285.9KB/s   00:13    
Fetching /home/tania/Popcorn/data2/H2-SIM-700026520-60.gwf to H2-SIM-700026520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700026520-60 100% 3733KB 155.5KB/s   00:24    
Fetching /home/tania/Popcorn/data2/H2-SIM-700026580-60.gwf to H2-SIM-700026580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700026580-60 100% 3715KB 742.9KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700026640-60.gwf to H2-SIM-700026640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700026640-60 100% 3708KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700026700-60.gwf to H2-SIM-700026700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700026700-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700026760-60.gwf to H2-SIM-700026760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700026760-60 100% 3706KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700026820-60.gwf to H2-SIM-700026820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700026820-60 100% 3713KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700026880-60.gwf to H2-SIM-700026880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700026880-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700026940-60.gwf to H2-SIM-700026940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700026940-60 100% 3704KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700027000-60.gwf to H2-SIM-700027000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700027000-60 100% 3705KB 463.2KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700027060-60.gwf to H2-SIM-700027060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700027060-60 100% 3708KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700027120-60.gwf to H2-SIM-700027120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700027120-60 100% 3701KB 740.1KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700027180-60.gwf to H2-SIM-700027180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700027180-60 100% 3717KB 743.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700027240-60.gwf to H2-SIM-700027240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700027240-60 100% 3710KB 463.7KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700027300-60.gwf to H2-SIM-700027300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700027300-60 100% 3716KB 119.9KB/s   00:31    
Fetching /home/tania/Popcorn/data2/H2-SIM-700027360-60.gwf to H2-SIM-700027360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700027360-60 100% 3713KB 742.6KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700027420-60.gwf to H2-SIM-700027420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700027420-60 100% 3715KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700027480-60.gwf to H2-SIM-700027480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700027480-60 100% 3712KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700027540-60.gwf to H2-SIM-700027540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700027540-60 100% 3713KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700027600-60.gwf to H2-SIM-700027600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700027600-60 100% 3709KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700027660-60.gwf to H2-SIM-700027660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700027660-60 100% 3712KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700027720-60.gwf to H2-SIM-700027720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700027720-60 100% 3711KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700027780-60.gwf to H2-SIM-700027780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700027780-60 100% 3708KB 370.8KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700027840-60.gwf to H2-SIM-700027840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700027840-60 100% 3710KB 142.7KB/s   00:26    
Fetching /home/tania/Popcorn/data2/H2-SIM-700027900-60.gwf to H2-SIM-700027900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700027900-60 100% 3704KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700027960-60.gwf to H2-SIM-700027960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700027960-60 100% 3716KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700028020-60.gwf to H2-SIM-700028020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700028020-60 100% 3699KB 739.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700028080-60.gwf to H2-SIM-700028080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700028080-60 100% 3713KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700028140-60.gwf to H2-SIM-700028140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700028140-60 100% 3701KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700028200-60.gwf to H2-SIM-700028200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700028200-60 100% 3706KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700028260-60.gwf to H2-SIM-700028260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700028260-60 100% 3703KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700028320-60.gwf to H2-SIM-700028320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700028320-60 100% 3708KB 741.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700028380-60.gwf to H2-SIM-700028380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700028380-60 100% 3707KB 926.7KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700028440-60.gwf to H2-SIM-700028440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700028440-60 100% 3712KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700028500-60.gwf to H2-SIM-700028500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700028500-60 100% 3717KB 285.9KB/s   00:13    
Fetching /home/tania/Popcorn/data2/H2-SIM-700028560-60.gwf to H2-SIM-700028560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700028560-60 100% 3722KB 124.1KB/s   00:30    
Fetching /home/tania/Popcorn/data2/H2-SIM-700028620-60.gwf to H2-SIM-700028620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700028620-60 100% 3712KB 742.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700028680-60.gwf to H2-SIM-700028680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700028680-60 100% 3713KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700028740-60.gwf to H2-SIM-700028740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700028740-60 100% 3703KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700028800-60.gwf to H2-SIM-700028800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700028800-60 100% 3715KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700028860-60.gwf to H2-SIM-700028860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700028860-60 100% 3712KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700028920-60.gwf to H2-SIM-700028920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700028920-60 100% 3709KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700028980-60.gwf to H2-SIM-700028980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700028980-60 100% 3715KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700029040-60.gwf to H2-SIM-700029040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700029040-60 100% 3710KB 337.3KB/s   00:11    
Fetching /home/tania/Popcorn/data2/H2-SIM-700029100-60.gwf to H2-SIM-700029100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700029100-60 100% 3722KB 128.3KB/s   00:29    
Fetching /home/tania/Popcorn/data2/H2-SIM-700029160-60.gwf to H2-SIM-700029160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700029160-60 100% 3722KB 465.2KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700029220-60.gwf to H2-SIM-700029220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700029220-60 100% 3711KB 927.7KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700029280-60.gwf to H2-SIM-700029280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700029280-60 100% 3713KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700029340-60.gwf to H2-SIM-700029340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700029340-60 100% 3707KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700029400-60.gwf to H2-SIM-700029400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700029400-60 100% 3711KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700029460-60.gwf to H2-SIM-700029460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700029460-60 100% 3709KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700029520-60.gwf to H2-SIM-700029520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700029520-60 100% 3705KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700029580-60.gwf to H2-SIM-700029580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700029580-60 100% 3703KB 740.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700029640-60.gwf to H2-SIM-700029640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700029640-60 100% 3714KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700029700-60.gwf to H2-SIM-700029700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700029700-60 100% 3702KB 740.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700029760-60.gwf to H2-SIM-700029760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700029760-60 100% 3716KB 530.8KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700029820-60.gwf to H2-SIM-700029820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700029820-60 100% 3703KB 529.1KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700029880-60.gwf to H2-SIM-700029880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700029880-60 100% 3710KB 618.3KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700029940-60.gwf to H2-SIM-700029940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700029940-60   0%    0     0.0KB/s   --:-- ETA
/home/tania/Popcorn/data2/H2-SIM-700029940-60 100% 3718KB 128.2KB/s   00:29    
Fetching /home/tania/Popcorn/data2/H2-SIM-700030000-60.gwf to H2-SIM-700030000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700030000-60 100% 3703KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700030060-60.gwf to H2-SIM-700030060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700030060-60 100% 3723KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700030120-60.gwf to H2-SIM-700030120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700030120-60 100% 3720KB 372.0KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700030180-60.gwf to H2-SIM-700030180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700030180-60 100% 3723KB 248.2KB/s   00:15    
Fetching /home/tania/Popcorn/data2/H2-SIM-700030240-60.gwf to H2-SIM-700030240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700030240-60 100% 3710KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700030300-60.gwf to H2-SIM-700030300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700030300-60  12%  448KB  13.1KB/s   04:08 ETA




/home/tania/Popcorn/data2/H2-SIM-700030300-60 100% 3710KB 218.2KB/s   00:17    
Fetching /home/tania/Popcorn/data2/H2-SIM-700030360-60.gwf to H2-SIM-700030360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700030360-60 100% 3707KB 411.9KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700030420-60.gwf to H2-SIM-700030420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700030420-60 100% 3708KB 741.6KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700030480-60.gwf to H2-SIM-700030480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700030480-60 100% 3705KB 463.2KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700030540-60.gwf to H2-SIM-700030540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700030540-60 100% 3709KB 148.4KB/s   00:25    
Fetching /home/tania/Popcorn/data2/H2-SIM-700030600-60.gwf to H2-SIM-700030600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700030600-60 100% 3698KB 924.6KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700030660-60.gwf to H2-SIM-700030660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700030660-60 100% 3725KB 745.0KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700030720-60.gwf to H2-SIM-700030720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700030720-60 100% 3716KB 743.2KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700030780-60.gwf to H2-SIM-700030780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700030780-60 100% 3707KB 264.8KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700030840-60.gwf to H2-SIM-700030840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700030840-60 100% 3721KB 177.2KB/s   00:21    
Fetching /home/tania/Popcorn/data2/H2-SIM-700030900-60.gwf to H2-SIM-700030900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700030900-60 100% 3711KB 927.7KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700030960-60.gwf to H2-SIM-700030960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700030960-60 100% 3700KB 739.9KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700031020-60.gwf to H2-SIM-700031020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700031020-60 100% 3714KB 530.6KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700031080-60.gwf to H2-SIM-700031080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700031080-60 100% 3708KB 154.5KB/s   00:24    
Fetching /home/tania/Popcorn/data2/H2-SIM-700031140-60.gwf to H2-SIM-700031140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700031140-60 100% 3718KB 619.6KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700031200-60.gwf to H2-SIM-700031200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700031200-60 100% 3702KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700031260-60.gwf to H2-SIM-700031260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700031260-60 100% 3704KB 463.0KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700031320-60.gwf to H2-SIM-700031320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700031320-60 100% 3706KB 161.1KB/s   00:23    
Fetching /home/tania/Popcorn/data2/H2-SIM-700031380-60.gwf to H2-SIM-700031380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700031380-60 100% 3710KB 927.5KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700031440-60.gwf to H2-SIM-700031440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700031440-60 100% 3723KB 620.5KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700031500-60.gwf to H2-SIM-700031500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700031500-60 100% 3710KB 412.3KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700031560-60.gwf to H2-SIM-700031560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700031560-60 100% 3714KB 218.5KB/s   00:17    
Fetching /home/tania/Popcorn/data2/H2-SIM-700031620-60.gwf to H2-SIM-700031620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700031620-60 100% 3709KB 412.1KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700031680-60.gwf to H2-SIM-700031680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700031680-60 100% 3710KB 412.2KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700031740-60.gwf to H2-SIM-700031740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700031740-60 100% 3712KB 142.8KB/s   00:26    
Fetching /home/tania/Popcorn/data2/H2-SIM-700031800-60.gwf to H2-SIM-700031800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700031800-60 100% 3704KB 617.4KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700031860-60.gwf to H2-SIM-700031860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700031860-60 100% 3699KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700031920-60.gwf to H2-SIM-700031920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700031920-60 100% 3714KB 464.2KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700031980-60.gwf to H2-SIM-700031980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700031980-60 100% 3712KB 218.3KB/s   00:17    
Fetching /home/tania/Popcorn/data2/H2-SIM-700032040-60.gwf to H2-SIM-700032040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700032040-60 100% 3718KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700032100-60.gwf to H2-SIM-700032100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700032100-60 100% 3707KB 218.1KB/s   00:17    
Fetching /home/tania/Popcorn/data2/H2-SIM-700032160-60.gwf to H2-SIM-700032160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700032160-60 100% 3713KB 337.5KB/s   00:11    
Fetching /home/tania/Popcorn/data2/H2-SIM-700032220-60.gwf to H2-SIM-700032220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700032220-60 100% 3690KB 922.4KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700032280-60.gwf to H2-SIM-700032280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700032280-60 100% 3710KB 742.1KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700032340-60.gwf to H2-SIM-700032340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700032340-60 100% 3707KB 264.8KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700032400-60.gwf to H2-SIM-700032400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700032400-60 100% 3722KB 169.2KB/s   00:22    
Fetching /home/tania/Popcorn/data2/H2-SIM-700032460-60.gwf to H2-SIM-700032460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700032460-60 100% 3697KB 739.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700032520-60.gwf to H2-SIM-700032520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700032520-60 100% 3707KB 529.6KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700032580-60.gwf to H2-SIM-700032580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700032580-60 100% 3715KB 176.9KB/s   00:21    
Fetching /home/tania/Popcorn/data2/H2-SIM-700032640-60.gwf to H2-SIM-700032640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700032640-60 100% 3705KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700032700-60.gwf to H2-SIM-700032700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700032700-60 100% 3709KB 741.9KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700032760-60.gwf to H2-SIM-700032760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700032760-60 100% 3707KB 741.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700032820-60.gwf to H2-SIM-700032820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700032820-60 100% 3710KB 247.3KB/s   00:15    
Fetching /home/tania/Popcorn/data2/H2-SIM-700032880-60.gwf to H2-SIM-700032880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700032880-60 100% 3698KB 137.0KB/s   00:27    
Fetching /home/tania/Popcorn/data2/H2-SIM-700032940-60.gwf to H2-SIM-700032940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700032940-60 100% 3715KB 619.2KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700033000-60.gwf to H2-SIM-700033000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700033000-60 100% 3702KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700033060-60.gwf to H2-SIM-700033060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700033060-60 100% 3702KB 740.4KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700033120-60.gwf to H2-SIM-700033120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700033120-60 100% 3708KB 264.8KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700033180-60.gwf to H2-SIM-700033180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700033180-60 100% 3721KB 161.8KB/s   00:23    
Fetching /home/tania/Popcorn/data2/H2-SIM-700033240-60.gwf to H2-SIM-700033240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700033240-60 100% 3701KB 616.8KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700033300-60.gwf to H2-SIM-700033300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700033300-60 100% 3708KB 741.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700033360-60.gwf to H2-SIM-700033360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700033360-60 100% 3716KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700033420-60.gwf to H2-SIM-700033420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700033420-60 100% 3716KB 177.0KB/s   00:21    
Fetching /home/tania/Popcorn/data2/H2-SIM-700033480-60.gwf to H2-SIM-700033480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700033480-60 100% 3717KB 206.5KB/s   00:18    
Fetching /home/tania/Popcorn/data2/H2-SIM-700033540-60.gwf to H2-SIM-700033540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700033540-60 100% 3709KB 741.7KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700033600-60.gwf to H2-SIM-700033600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700033600-60 100% 3709KB 618.2KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700033660-60.gwf to H2-SIM-700033660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700033660-60 100% 3710KB 742.0KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700033720-60.gwf to H2-SIM-700033720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700033720-60 100% 3707KB 176.5KB/s   00:21    
Fetching /home/tania/Popcorn/data2/H2-SIM-700033780-60.gwf to H2-SIM-700033780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700033780-60 100% 3714KB 371.4KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700033840-60.gwf to H2-SIM-700033840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700033840-60 100% 3695KB 461.9KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700033900-60.gwf to H2-SIM-700033900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700033900-60 100% 3707KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700033960-60.gwf to H2-SIM-700033960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700033960-60 100% 3720KB 248.0KB/s   00:15    
Fetching /home/tania/Popcorn/data2/H2-SIM-700034020-60.gwf to H2-SIM-700034020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700034020-60 100% 3723KB 143.2KB/s   00:26    
Fetching /home/tania/Popcorn/data2/H2-SIM-700034080-60.gwf to H2-SIM-700034080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700034080-60 100% 3718KB 929.6KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700034140-60.gwf to H2-SIM-700034140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700034140-60 100% 3711KB 463.9KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700034200-60.gwf to H2-SIM-700034200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700034200-60 100% 3704KB 308.7KB/s   00:12    
Fetching /home/tania/Popcorn/data2/H2-SIM-700034260-60.gwf to H2-SIM-700034260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700034260-60 100% 3708KB 142.6KB/s   00:26    
Fetching /home/tania/Popcorn/data2/H2-SIM-700034320-60.gwf to H2-SIM-700034320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700034320-60 100% 3703KB 617.2KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700034380-60.gwf to H2-SIM-700034380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700034380-60 100% 3699KB 616.5KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700034440-60.gwf to H2-SIM-700034440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700034440-60 100% 3705KB 741.0KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700034500-60.gwf to H2-SIM-700034500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700034500-60   0%    0     0.0KB/s   --:-- ETA

/home/tania/Popcorn/data2/H2-SIM-700034500-60 100% 3705KB 142.5KB/s   00:26    
Fetching /home/tania/Popcorn/data2/H2-SIM-700034560-60.gwf to H2-SIM-700034560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700034560-60 100% 3708KB 370.8KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700034620-60.gwf to H2-SIM-700034620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700034620-60 100% 3711KB 742.2KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700034680-60.gwf to H2-SIM-700034680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700034680-60 100% 3716KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700034740-60.gwf to H2-SIM-700034740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700034740-60 100% 3718KB 413.1KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700034800-60.gwf to H2-SIM-700034800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700034800-60 100% 3714KB 137.6KB/s   00:27    
Fetching /home/tania/Popcorn/data2/H2-SIM-700034860-60.gwf to H2-SIM-700034860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700034860-60 100% 3705KB 411.7KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700034920-60.gwf to H2-SIM-700034920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700034920-60 100% 3714KB 742.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700034980-60.gwf to H2-SIM-700034980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700034980-60 100% 3711KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700035040-60.gwf to H2-SIM-700035040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700035040-60 100% 3706KB 247.1KB/s   00:15    
Fetching /home/tania/Popcorn/data2/H2-SIM-700035100-60.gwf to H2-SIM-700035100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700035100-60 100% 3706KB 168.4KB/s   00:22    
Fetching /home/tania/Popcorn/data2/H2-SIM-700035160-60.gwf to H2-SIM-700035160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700035160-60 100% 3713KB 371.4KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700035220-60.gwf to H2-SIM-700035220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700035220-60 100% 3713KB 742.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700035280-60.gwf to H2-SIM-700035280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700035280-60 100% 3708KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700035340-60.gwf to H2-SIM-700035340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700035340-60 100% 3710KB 371.0KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700035400-60.gwf to H2-SIM-700035400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700035400-60 100% 3697KB 142.2KB/s   00:26    
Fetching /home/tania/Popcorn/data2/H2-SIM-700035460-60.gwf to H2-SIM-700035460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700035460-60 100% 3714KB 530.6KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700035520-60.gwf to H2-SIM-700035520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700035520-60 100% 3725KB 744.9KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700035580-60.gwf to H2-SIM-700035580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700035580-60 100% 3708KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700035640-60.gwf to H2-SIM-700035640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700035640-60 100% 3699KB 205.5KB/s   00:18    
Fetching /home/tania/Popcorn/data2/H2-SIM-700035700-60.gwf to H2-SIM-700035700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700035700-60 100% 3708KB 176.6KB/s   00:21    
Fetching /home/tania/Popcorn/data2/H2-SIM-700035760-60.gwf to H2-SIM-700035760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700035760-60 100% 3717KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700035820-60.gwf to H2-SIM-700035820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700035820-60 100% 3713KB 618.8KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700035880-60.gwf to H2-SIM-700035880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700035880-60 100% 3713KB 928.2KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700035940-60.gwf to H2-SIM-700035940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700035940-60 100% 3716KB 218.6KB/s   00:17    
Fetching /home/tania/Popcorn/data2/H2-SIM-700036000-60.gwf to H2-SIM-700036000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700036000-60 100% 3712KB 161.4KB/s   00:23    
Fetching /home/tania/Popcorn/data2/H2-SIM-700036060-60.gwf to H2-SIM-700036060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700036060-60 100% 3711KB 742.2KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700036120-60.gwf to H2-SIM-700036120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700036120-60 100% 3706KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700036180-60.gwf to H2-SIM-700036180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700036180-60 100% 3704KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700036240-60.gwf to H2-SIM-700036240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700036240-60 100% 3715KB 195.5KB/s   00:19    
Fetching /home/tania/Popcorn/data2/H2-SIM-700036300-60.gwf to H2-SIM-700036300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700036300-60 100% 3710KB 218.3KB/s   00:17    
Fetching /home/tania/Popcorn/data2/H2-SIM-700036360-60.gwf to H2-SIM-700036360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700036360-60 100% 3727KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700036420-60.gwf to H2-SIM-700036420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700036420-60 100% 3718KB 743.6KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700036480-60.gwf to H2-SIM-700036480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700036480-60 100% 3699KB 369.9KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700036540-60.gwf to H2-SIM-700036540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700036540-60 100% 3716KB 177.0KB/s   00:21    
Fetching /home/tania/Popcorn/data2/H2-SIM-700036600-60.gwf to H2-SIM-700036600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700036600-60 100% 3714KB 530.6KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700036660-60.gwf to H2-SIM-700036660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700036660-60 100% 3715KB 206.4KB/s   00:18    
Fetching /home/tania/Popcorn/data2/H2-SIM-700036720-60.gwf to H2-SIM-700036720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700036720-60 100% 3713KB 185.7KB/s   00:20    
Fetching /home/tania/Popcorn/data2/H2-SIM-700036780-60.gwf to H2-SIM-700036780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700036780-60 100% 3708KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700036840-60.gwf to H2-SIM-700036840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700036840-60 100% 3718KB 929.4KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700036900-60.gwf to H2-SIM-700036900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700036900-60 100% 3704KB 246.9KB/s   00:15    
Fetching /home/tania/Popcorn/data2/H2-SIM-700036960-60.gwf to H2-SIM-700036960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700036960-60 100% 3709KB 161.3KB/s   00:23    
Fetching /home/tania/Popcorn/data2/H2-SIM-700037020-60.gwf to H2-SIM-700037020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700037020-60 100% 3711KB 742.2KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700037080-60.gwf to H2-SIM-700037080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700037080-60 100% 3689KB 737.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700037140-60.gwf to H2-SIM-700037140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700037140-60 100% 3702KB 462.8KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700037200-60.gwf to H2-SIM-700037200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700037200-60 100% 3713KB 128.0KB/s   00:29    
Fetching /home/tania/Popcorn/data2/H2-SIM-700037260-60.gwf to H2-SIM-700037260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700037260-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700037260-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700037320-60.gwf to H2-SIM-700037320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700037320-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700037320-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700037380-60.gwf to H2-SIM-700037380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700037380-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700037380-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700037440-60.gwf to H2-SIM-700037440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700037440-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700037440-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700037500-60.gwf to H2-SIM-700037500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700037500-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700037500-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700037560-60.gwf to H2-SIM-700037560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700037560-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700037560-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700037620-60.gwf to H2-SIM-700037620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700037620-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700037620-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700037680-60.gwf to H2-SIM-700037680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700037680-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700037680-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700037740-60.gwf to H2-SIM-700037740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700037740-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700037740-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700037800-60.gwf to H2-SIM-700037800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700037800-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700037800-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700037860-60.gwf to H2-SIM-700037860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700037860-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700037860-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700037920-60.gwf to H2-SIM-700037920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700037920-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700037920-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700037980-60.gwf to H2-SIM-700037980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700037980-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700037980-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700038040-60.gwf to H2-SIM-700038040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700038040-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700038040-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700038100-60.gwf to H2-SIM-700038100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700038100-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700038100-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700038160-60.gwf to H2-SIM-700038160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700038160-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700038160-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700038220-60.gwf to H2-SIM-700038220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700038220-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700038220-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700038280-60.gwf to H2-SIM-700038280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700038280-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700038280-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700038340-60.gwf to H2-SIM-700038340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700038340-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700038340-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700038400-60.gwf to H2-SIM-700038400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700038400-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700038400-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700038460-60.gwf to H2-SIM-700038460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700038460-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700038460-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700038520-60.gwf to H2-SIM-700038520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700038520-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700038520-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700038580-60.gwf to H2-SIM-700038580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700038580-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700038580-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700038640-60.gwf to H2-SIM-700038640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700038640-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700038640-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700038700-60.gwf to H2-SIM-700038700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700038700-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700038700-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700038760-60.gwf to H2-SIM-700038760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700038760-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700038760-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700038820-60.gwf to H2-SIM-700038820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700038820-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700038820-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700038880-60.gwf to H2-SIM-700038880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700038880-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700038880-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700038940-60.gwf to H2-SIM-700038940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700038940-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700038940-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700039000-60.gwf to H2-SIM-700039000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700039000-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700039000-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700039060-60.gwf to H2-SIM-700039060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700039060-60 100% 3723KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700039120-60.gwf to H2-SIM-700039120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700039120-60 100% 3707KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700039180-60.gwf to H2-SIM-700039180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700039180-60 100% 3702KB 411.4KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700039240-60.gwf to H2-SIM-700039240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700039240-60 100% 3702KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700039300-60.gwf to H2-SIM-700039300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700039300-60 100% 3703KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700039360-60.gwf to H2-SIM-700039360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700039360-60 100% 3701KB 462.6KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700039420-60.gwf to H2-SIM-700039420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700039420-60 100% 3706KB 142.5KB/s   00:26    
Fetching /home/tania/Popcorn/data2/H2-SIM-700039480-60.gwf to H2-SIM-700039480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700039480-60 100% 3706KB 741.1KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700039540-60.gwf to H2-SIM-700039540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700039540-60 100% 3712KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700039600-60.gwf to H2-SIM-700039600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700039600-60 100% 3723KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700039660-60.gwf to H2-SIM-700039660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700039660-60 100% 3711KB 247.4KB/s   00:15    
Fetching /home/tania/Popcorn/data2/H2-SIM-700039720-60.gwf to H2-SIM-700039720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700039720-60 100% 3712KB 168.7KB/s   00:22    
Fetching /home/tania/Popcorn/data2/H2-SIM-700039780-60.gwf to H2-SIM-700039780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700039780-60 100% 3706KB 529.4KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700039840-60.gwf to H2-SIM-700039840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700039840-60 100% 3716KB 743.2KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700039900-60.gwf to H2-SIM-700039900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700039900-60 100% 3705KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700039960-60.gwf to H2-SIM-700039960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700039960-60 100% 3706KB 617.6KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700040020-60.gwf to H2-SIM-700040020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700040020-60 100% 3710KB 119.7KB/s   00:31    
Fetching /home/tania/Popcorn/data2/H2-SIM-700040080-60.gwf to H2-SIM-700040080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700040080-60 100% 3720KB 286.1KB/s   00:13    
Fetching /home/tania/Popcorn/data2/H2-SIM-700040140-60.gwf to H2-SIM-700040140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700040140-60 100% 3703KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700040200-60.gwf to H2-SIM-700040200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700040200-60 100% 3702KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700040260-60.gwf to H2-SIM-700040260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700040260-60 100% 3715KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700040320-60.gwf to H2-SIM-700040320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700040320-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700040380-60.gwf to H2-SIM-700040380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700040380-60 100% 3707KB 617.8KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700040440-60.gwf to H2-SIM-700040440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700040440-60 100% 3700KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700040500-60.gwf to H2-SIM-700040500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700040500-60 100% 3706KB 926.6KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700040560-60.gwf to H2-SIM-700040560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700040560-60 100% 3712KB 412.5KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700040620-60.gwf to H2-SIM-700040620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700040620-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700040620-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700040680-60.gwf to H2-SIM-700040680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700040680-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700040680-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700040740-60.gwf to H2-SIM-700040740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700040740-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700040740-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700040800-60.gwf to H2-SIM-700040800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700040800-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700040800-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700040860-60.gwf to H2-SIM-700040860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700040860-60   0%   32KB  32.0KB/s   01:55 ETA
Couldn't write to "H2-SIM-700040860-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700040920-60.gwf to H2-SIM-700040920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700040920-60   0%   32KB  32.0KB/s   01:54 ETA
Couldn't write to "H2-SIM-700040920-60.gwf": Input/output error
Fetching /home/tania/Popcorn/data2/H2-SIM-700040980-60.gwf to H2-SIM-700040980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700040980-60 100% 3696KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700041040-60.gwf to H2-SIM-700041040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700041040-60 100% 3727KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700041100-60.gwf to H2-SIM-700041100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700041100-60 100% 3703KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700041160-60.gwf to H2-SIM-700041160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700041160-60 100% 3720KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700041220-60.gwf to H2-SIM-700041220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700041220-60 100% 3700KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700041280-60.gwf to H2-SIM-700041280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700041280-60 100% 3712KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700041340-60.gwf to H2-SIM-700041340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700041340-60 100% 3706KB 926.5KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700041400-60.gwf to H2-SIM-700041400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700041400-60 100% 3707KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700041460-60.gwf to H2-SIM-700041460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700041460-60 100% 3710KB 927.4KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700041520-60.gwf to H2-SIM-700041520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700041520-60 100% 3706KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700041580-60.gwf to H2-SIM-700041580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700041580-60 100% 3699KB 739.7KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700041640-60.gwf to H2-SIM-700041640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700041640-60 100% 3703KB 740.7KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700041700-60.gwf to H2-SIM-700041700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700041700-60 100% 3718KB 929.4KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700041760-60.gwf to H2-SIM-700041760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700041760-60 100% 3703KB 411.4KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700041820-60.gwf to H2-SIM-700041820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700041820-60 100% 3707KB 161.2KB/s   00:23    
Fetching /home/tania/Popcorn/data2/H2-SIM-700041880-60.gwf to H2-SIM-700041880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700041880-60 100% 3707KB 741.4KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700041940-60.gwf to H2-SIM-700041940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700041940-60 100% 3714KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700042000-60.gwf to H2-SIM-700042000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700042000-60 100% 3707KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700042060-60.gwf to H2-SIM-700042060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700042060-60 100% 3700KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700042120-60.gwf to H2-SIM-700042120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700042120-60 100% 3704KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700042180-60.gwf to H2-SIM-700042180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700042180-60 100% 3714KB 742.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700042240-60.gwf to H2-SIM-700042240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700042240-60 100% 3705KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700042300-60.gwf to H2-SIM-700042300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700042300-60 100% 3702KB 370.2KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700042360-60.gwf to H2-SIM-700042360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700042360-60 100% 3714KB 168.8KB/s   00:22    
Fetching /home/tania/Popcorn/data2/H2-SIM-700042420-60.gwf to H2-SIM-700042420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700042420-60 100% 3708KB 529.7KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700042480-60.gwf to H2-SIM-700042480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700042480-60 100% 3717KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700042540-60.gwf to H2-SIM-700042540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700042540-60 100% 3713KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700042600-60.gwf to H2-SIM-700042600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700042600-60 100% 3703KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700042660-60.gwf to H2-SIM-700042660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700042660-60 100% 3713KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700042720-60.gwf to H2-SIM-700042720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700042720-60 100% 3707KB 741.4KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700042780-60.gwf to H2-SIM-700042780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700042780-60 100% 3717KB 619.5KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700042840-60.gwf to H2-SIM-700042840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700042840-60 100% 3715KB 743.0KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700042900-60.gwf to H2-SIM-700042900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700042900-60 100% 3712KB 927.9KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700042960-60.gwf to H2-SIM-700042960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700042960-60 100% 3703KB 370.3KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700043020-60.gwf to H2-SIM-700043020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700043020-60 100% 3710KB 123.7KB/s   00:30    
Fetching /home/tania/Popcorn/data2/H2-SIM-700043080-60.gwf to H2-SIM-700043080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700043080-60 100% 3713KB 928.1KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700043140-60.gwf to H2-SIM-700043140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700043140-60 100% 3711KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700043200-60.gwf to H2-SIM-700043200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700043200-60 100% 3707KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700043260-60.gwf to H2-SIM-700043260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700043260-60 100% 3718KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700043320-60.gwf to H2-SIM-700043320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700043320-60 100% 3712KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700043380-60.gwf to H2-SIM-700043380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700043380-60 100% 3715KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700043440-60.gwf to H2-SIM-700043440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700043440-60 100% 3704KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700043500-60.gwf to H2-SIM-700043500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700043500-60 100% 3706KB 926.6KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700043560-60.gwf to H2-SIM-700043560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700043560-60 100% 3701KB 137.1KB/s   00:27    
Fetching /home/tania/Popcorn/data2/H2-SIM-700043620-60.gwf to H2-SIM-700043620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700043620-60 100% 3712KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700043680-60.gwf to H2-SIM-700043680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700043680-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700043740-60.gwf to H2-SIM-700043740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700043740-60 100% 3695KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700043800-60.gwf to H2-SIM-700043800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700043800-60 100% 3709KB 618.2KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700043860-60.gwf to H2-SIM-700043860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700043860-60 100% 3718KB 929.5KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700043920-60.gwf to H2-SIM-700043920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700043920-60 100% 3711KB 742.1KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700043980-60.gwf to H2-SIM-700043980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700043980-60 100% 3708KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700044040-60.gwf to H2-SIM-700044040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700044040-60 100% 3709KB 927.2KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700044100-60.gwf to H2-SIM-700044100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700044100-60 100% 3714KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700044160-60.gwf to H2-SIM-700044160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700044160-60 100% 3707KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700044220-60.gwf to H2-SIM-700044220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700044220-60 100% 3712KB 337.5KB/s   00:11    
Fetching /home/tania/Popcorn/data2/H2-SIM-700044280-60.gwf to H2-SIM-700044280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700044280-60 100% 3710KB 112.4KB/s   00:33    
Fetching /home/tania/Popcorn/data2/H2-SIM-700044340-60.gwf to H2-SIM-700044340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700044340-60 100% 3708KB 309.0KB/s   00:12    
Fetching /home/tania/Popcorn/data2/H2-SIM-700044400-60.gwf to H2-SIM-700044400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700044400-60 100% 3714KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700044460-60.gwf to H2-SIM-700044460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700044460-60 100% 3720KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700044520-60.gwf to H2-SIM-700044520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700044520-60 100% 3706KB 741.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700044580-60.gwf to H2-SIM-700044580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700044580-60 100% 3716KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700044640-60.gwf to H2-SIM-700044640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700044640-60 100% 3722KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700044700-60.gwf to H2-SIM-700044700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700044700-60 100% 3703KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700044760-60.gwf to H2-SIM-700044760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700044760-60 100% 3701KB 925.4KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700044820-60.gwf to H2-SIM-700044820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700044820-60 100% 3703KB 925.7KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700044880-60.gwf to H2-SIM-700044880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700044880-60 100% 3709KB 412.2KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700044940-60.gwf to H2-SIM-700044940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700044940-60 100% 3703KB 127.7KB/s   00:29    
Fetching /home/tania/Popcorn/data2/H2-SIM-700045000-60.gwf to H2-SIM-700045000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700045000-60 100% 3709KB 741.9KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700045060-60.gwf to H2-SIM-700045060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700045060-60 100% 3717KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700045120-60.gwf to H2-SIM-700045120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700045120-60 100% 3700KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700045180-60.gwf to H2-SIM-700045180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700045180-60 100% 3719KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700045240-60.gwf to H2-SIM-700045240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700045240-60 100% 3723KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700045300-60.gwf to H2-SIM-700045300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700045300-60 100% 3712KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700045360-60.gwf to H2-SIM-700045360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700045360-60 100% 3696KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700045420-60.gwf to H2-SIM-700045420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700045420-60 100% 3706KB 926.5KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700045480-60.gwf to H2-SIM-700045480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700045480-60 100% 3707KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700045540-60.gwf to H2-SIM-700045540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700045540-60 100% 3715KB 619.2KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700045600-60.gwf to H2-SIM-700045600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700045600-60 100% 3713KB 742.7KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700045660-60.gwf to H2-SIM-700045660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700045660-60 100% 3705KB 463.2KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700045720-60.gwf to H2-SIM-700045720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700045720-60 100% 3710KB 123.7KB/s   00:30    
Fetching /home/tania/Popcorn/data2/H2-SIM-700045780-60.gwf to H2-SIM-700045780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700045780-60 100% 3713KB 742.6KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700045840-60.gwf to H2-SIM-700045840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700045840-60 100% 3714KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700045900-60.gwf to H2-SIM-700045900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700045900-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700045960-60.gwf to H2-SIM-700045960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700045960-60 100% 3705KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700046020-60.gwf to H2-SIM-700046020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700046020-60 100% 3711KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700046080-60.gwf to H2-SIM-700046080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700046080-60 100% 3705KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700046140-60.gwf to H2-SIM-700046140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700046140-60 100% 3700KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700046200-60.gwf to H2-SIM-700046200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700046200-60 100% 3714KB 285.7KB/s   00:13    
Fetching /home/tania/Popcorn/data2/H2-SIM-700046260-60.gwf to H2-SIM-700046260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700046260-60 100% 3719KB 155.0KB/s   00:24    
Fetching /home/tania/Popcorn/data2/H2-SIM-700046320-60.gwf to H2-SIM-700046320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700046320-60 100% 3708KB 412.0KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700046380-60.gwf to H2-SIM-700046380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700046380-60 100% 3718KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700046440-60.gwf to H2-SIM-700046440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700046440-60 100% 3712KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700046500-60.gwf to H2-SIM-700046500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700046500-60 100% 3717KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700046560-60.gwf to H2-SIM-700046560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700046560-60 100% 3726KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700046620-60.gwf to H2-SIM-700046620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700046620-60 100% 3704KB 740.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700046680-60.gwf to H2-SIM-700046680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700046680-60 100% 3718KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700046740-60.gwf to H2-SIM-700046740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700046740-60 100% 3724KB 744.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700046800-60.gwf to H2-SIM-700046800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700046800-60 100% 3714KB 928.6KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700046860-60.gwf to H2-SIM-700046860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700046860-60 100% 3708KB 741.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700046920-60.gwf to H2-SIM-700046920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700046920-60 100% 3715KB 265.4KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700046980-60.gwf to H2-SIM-700046980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700046980-60 100% 3723KB 143.2KB/s   00:26    
Fetching /home/tania/Popcorn/data2/H2-SIM-700047040-60.gwf to H2-SIM-700047040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700047040-60 100% 3701KB 740.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700047100-60.gwf to H2-SIM-700047100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700047100-60 100% 3699KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700047160-60.gwf to H2-SIM-700047160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700047160-60 100% 3710KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700047220-60.gwf to H2-SIM-700047220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700047220-60 100% 3710KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700047280-60.gwf to H2-SIM-700047280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700047280-60 100% 3712KB 928.1KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700047340-60.gwf to H2-SIM-700047340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700047340-60 100% 3704KB 185.2KB/s   00:20    
Fetching /home/tania/Popcorn/data2/H2-SIM-700047400-60.gwf to H2-SIM-700047400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700047400-60 100% 3715KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700047460-60.gwf to H2-SIM-700047460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700047460-60 100% 3720KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700047520-60.gwf to H2-SIM-700047520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700047520-60 100% 3701KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700047580-60.gwf to H2-SIM-700047580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700047580-60 100% 3707KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700047640-60.gwf to H2-SIM-700047640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700047640-60 100% 3716KB 743.1KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700047700-60.gwf to H2-SIM-700047700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700047700-60 100% 3712KB 265.1KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700047760-60.gwf to H2-SIM-700047760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700047760-60 100% 3716KB 218.6KB/s   00:17    
Fetching /home/tania/Popcorn/data2/H2-SIM-700047820-60.gwf to H2-SIM-700047820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700047820-60 100% 3705KB 741.0KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700047880-60.gwf to H2-SIM-700047880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700047880-60 100% 3706KB 926.5KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700047940-60.gwf to H2-SIM-700047940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700047940-60 100% 3706KB 741.2KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700048000-60.gwf to H2-SIM-700048000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700048000-60 100% 3701KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700048060-60.gwf to H2-SIM-700048060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700048060-60 100% 3702KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700048120-60.gwf to H2-SIM-700048120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700048120-60 100% 3711KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700048180-60.gwf to H2-SIM-700048180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700048180-60 100% 3713KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700048240-60.gwf to H2-SIM-700048240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700048240-60 100% 3709KB 206.0KB/s   00:18    
Fetching /home/tania/Popcorn/data2/H2-SIM-700048300-60.gwf to H2-SIM-700048300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700048300-60 100% 3716KB 743.1KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700048360-60.gwf to H2-SIM-700048360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700048360-60 100% 3698KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700048420-60.gwf to H2-SIM-700048420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700048420-60 100% 3708KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700048480-60.gwf to H2-SIM-700048480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700048480-60 100% 3716KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700048540-60.gwf to H2-SIM-700048540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700048540-60 100% 3720KB 218.8KB/s   00:17    
Fetching /home/tania/Popcorn/data2/H2-SIM-700048600-60.gwf to H2-SIM-700048600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700048600-60 100% 3713KB 195.4KB/s   00:19    
Fetching /home/tania/Popcorn/data2/H2-SIM-700048660-60.gwf to H2-SIM-700048660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700048660-60 100% 3717KB 743.4KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700048720-60.gwf to H2-SIM-700048720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700048720-60 100% 3714KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700048780-60.gwf to H2-SIM-700048780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700048780-60 100% 3699KB 924.6KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700048840-60.gwf to H2-SIM-700048840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700048840-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700048900-60.gwf to H2-SIM-700048900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700048900-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700048960-60.gwf to H2-SIM-700048960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700048960-60 100% 3716KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700049020-60.gwf to H2-SIM-700049020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700049020-60 100% 3715KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700049080-60.gwf to H2-SIM-700049080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700049080-60 100% 3723KB 744.6KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700049140-60.gwf to H2-SIM-700049140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700049140-60 100% 3709KB 927.2KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700049200-60.gwf to H2-SIM-700049200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700049200-60 100% 3708KB 618.0KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700049260-60.gwf to H2-SIM-700049260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700049260-60 100% 3711KB 265.1KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700049320-60.gwf to H2-SIM-700049320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700049320-60 100% 3702KB 142.4KB/s   00:26    
Fetching /home/tania/Popcorn/data2/H2-SIM-700049380-60.gwf to H2-SIM-700049380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700049380-60 100% 3707KB 741.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700049440-60.gwf to H2-SIM-700049440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700049440-60 100% 3708KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700049500-60.gwf to H2-SIM-700049500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700049500-60 100% 3695KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700049560-60.gwf to H2-SIM-700049560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700049560-60 100% 3698KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700049620-60.gwf to H2-SIM-700049620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700049620-60 100% 3718KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700049680-60.gwf to H2-SIM-700049680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700049680-60 100% 3710KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700049740-60.gwf to H2-SIM-700049740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700049740-60 100% 3712KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700049800-60.gwf to H2-SIM-700049800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700049800-60 100% 3714KB 619.1KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700049860-60.gwf to H2-SIM-700049860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700049860-60 100% 3709KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700049920-60.gwf to H2-SIM-700049920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700049920-60 100% 3715KB 928.8KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700049980-60.gwf to H2-SIM-700049980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700049980-60 100% 3717KB 371.7KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700050040-60.gwf to H2-SIM-700050040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700050040-60 100% 3722KB 248.2KB/s   00:15    
Fetching /home/tania/Popcorn/data2/H2-SIM-700050100-60.gwf to H2-SIM-700050100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700050100-60 100% 3711KB 142.7KB/s   00:26    
Fetching /home/tania/Popcorn/data2/H2-SIM-700050160-60.gwf to H2-SIM-700050160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700050160-60 100% 3710KB 742.0KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700050220-60.gwf to H2-SIM-700050220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700050220-60 100% 3713KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700050280-60.gwf to H2-SIM-700050280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700050280-60 100% 3713KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700050340-60.gwf to H2-SIM-700050340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700050340-60 100% 3707KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700050400-60.gwf to H2-SIM-700050400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700050400-60 100% 3714KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700050460-60.gwf to H2-SIM-700050460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700050460-60 100% 3710KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700050520-60.gwf to H2-SIM-700050520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700050520-60 100% 3704KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700050580-60.gwf to H2-SIM-700050580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700050580-60 100% 3706KB 926.5KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700050640-60.gwf to H2-SIM-700050640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700050640-60 100% 3708KB 463.5KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700050700-60.gwf to H2-SIM-700050700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700050700-60 100% 3719KB 743.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700050760-60.gwf to H2-SIM-700050760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700050760-60 100% 3712KB 928.0KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700050820-60.gwf to H2-SIM-700050820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700050820-60 100% 3697KB 231.0KB/s   00:16    
Fetching /home/tania/Popcorn/data2/H2-SIM-700050880-60.gwf to H2-SIM-700050880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700050880-60 100% 3708KB 148.3KB/s   00:25    
Fetching /home/tania/Popcorn/data2/H2-SIM-700050940-60.gwf to H2-SIM-700050940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700050940-60 100% 3714KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700051000-60.gwf to H2-SIM-700051000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700051000-60 100% 3705KB 926.1KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700051060-60.gwf to H2-SIM-700051060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700051060-60 100% 3705KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700051120-60.gwf to H2-SIM-700051120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700051120-60 100% 3705KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700051180-60.gwf to H2-SIM-700051180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700051180-60 100% 3710KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700051240-60.gwf to H2-SIM-700051240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700051240-60 100% 3720KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700051300-60.gwf to H2-SIM-700051300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700051300-60 100% 3718KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700051360-60.gwf to H2-SIM-700051360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700051360-60 100% 3705KB 741.0KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700051420-60.gwf to H2-SIM-700051420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700051420-60 100% 3715KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700051480-60.gwf to H2-SIM-700051480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700051480-60 100% 3713KB 618.8KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700051540-60.gwf to H2-SIM-700051540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700051540-60 100% 3710KB 412.2KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700051600-60.gwf to H2-SIM-700051600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700051600-60 100% 3703KB 132.3KB/s   00:28    
Fetching /home/tania/Popcorn/data2/H2-SIM-700051660-60.gwf to H2-SIM-700051660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700051660-60 100% 3716KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700051720-60.gwf to H2-SIM-700051720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700051720-60 100% 3694KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700051780-60.gwf to H2-SIM-700051780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700051780-60 100% 3701KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700051840-60.gwf to H2-SIM-700051840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700051840-60 100% 3714KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700051900-60.gwf to H2-SIM-700051900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700051900-60 100% 3718KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700051960-60.gwf to H2-SIM-700051960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700051960-60 100% 3708KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700052020-60.gwf to H2-SIM-700052020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700052020-60 100% 3713KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700052080-60.gwf to H2-SIM-700052080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700052080-60 100% 3710KB 742.0KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700052140-60.gwf to H2-SIM-700052140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700052140-60 100% 3706KB 926.6KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700052200-60.gwf to H2-SIM-700052200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700052200-60 100% 3710KB 742.0KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700052260-60.gwf to H2-SIM-700052260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700052260-60 100% 3726KB 745.2KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700052320-60.gwf to H2-SIM-700052320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700052320-60 100% 3700KB 411.1KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700052380-60.gwf to H2-SIM-700052380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700052380-60 100% 3709KB 127.9KB/s   00:29    
Fetching /home/tania/Popcorn/data2/H2-SIM-700052440-60.gwf to H2-SIM-700052440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700052440-60 100% 3694KB 738.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700052500-60.gwf to H2-SIM-700052500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700052500-60 100% 3704KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700052560-60.gwf to H2-SIM-700052560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700052560-60 100% 3709KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700052620-60.gwf to H2-SIM-700052620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700052620-60 100% 3713KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700052680-60.gwf to H2-SIM-700052680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700052680-60 100% 3709KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700052740-60.gwf to H2-SIM-700052740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700052740-60 100% 3711KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700052800-60.gwf to H2-SIM-700052800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700052800-60 100% 3700KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700052860-60.gwf to H2-SIM-700052860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700052860-60 100% 3699KB 336.3KB/s   00:11    
Fetching /home/tania/Popcorn/data2/H2-SIM-700052920-60.gwf to H2-SIM-700052920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700052920-60 100% 3714KB 337.6KB/s   00:11    
Fetching /home/tania/Popcorn/data2/H2-SIM-700052980-60.gwf to H2-SIM-700052980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700052980-60 100% 3716KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700053040-60.gwf to H2-SIM-700053040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700053040-60 100% 3716KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700053100-60.gwf to H2-SIM-700053100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700053100-60 100% 3703KB 740.6KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700053160-60.gwf to H2-SIM-700053160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700053160-60 100% 3707KB 741.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700053220-60.gwf to H2-SIM-700053220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700053220-60 100% 3704KB 740.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700053280-60.gwf to H2-SIM-700053280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700053280-60 100% 3709KB 412.1KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700053340-60.gwf to H2-SIM-700053340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700053340-60 100% 3703KB 154.3KB/s   00:24    
Fetching /home/tania/Popcorn/data2/H2-SIM-700053400-60.gwf to H2-SIM-700053400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700053400-60 100% 3700KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700053460-60.gwf to H2-SIM-700053460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700053460-60 100% 3720KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700053520-60.gwf to H2-SIM-700053520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700053520-60 100% 3697KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700053580-60.gwf to H2-SIM-700053580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700053580-60 100% 3710KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700053640-60.gwf to H2-SIM-700053640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700053640-60 100% 3696KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700053700-60.gwf to H2-SIM-700053700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700053700-60 100% 3715KB 743.0KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700053760-60.gwf to H2-SIM-700053760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700053760-60 100% 3714KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700053820-60.gwf to H2-SIM-700053820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700053820-60 100% 3714KB 928.5KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700053880-60.gwf to H2-SIM-700053880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700053880-60 100% 3697KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700053940-60.gwf to H2-SIM-700053940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700053940-60 100% 3711KB 927.8KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700054000-60.gwf to H2-SIM-700054000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700054000-60 100% 3704KB 168.4KB/s   00:22    
Fetching /home/tania/Popcorn/data2/H2-SIM-700054060-60.gwf to H2-SIM-700054060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700054060-60 100% 3715KB 154.8KB/s   00:24    
Fetching /home/tania/Popcorn/data2/H2-SIM-700054120-60.gwf to H2-SIM-700054120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700054120-60 100% 3707KB 741.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700054180-60.gwf to H2-SIM-700054180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700054180-60 100% 3722KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700054240-60.gwf to H2-SIM-700054240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700054240-60 100% 3703KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700054300-60.gwf to H2-SIM-700054300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700054300-60 100% 3705KB 926.3KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700054360-60.gwf to H2-SIM-700054360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700054360-60 100% 3711KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700054420-60.gwf to H2-SIM-700054420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700054420-60 100% 3714KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700054480-60.gwf to H2-SIM-700054480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700054480-60 100% 3707KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700054540-60.gwf to H2-SIM-700054540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700054540-60 100% 3711KB 742.2KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700054600-60.gwf to H2-SIM-700054600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700054600-60 100% 3718KB 929.4KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700054660-60.gwf to H2-SIM-700054660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700054660-60 100% 3709KB 741.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700054720-60.gwf to H2-SIM-700054720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700054720-60 100% 3704KB 740.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700054780-60.gwf to H2-SIM-700054780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700054780-60 100% 3706KB 741.1KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700054840-60.gwf to H2-SIM-700054840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700054840-60 100% 3717KB 464.6KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700054900-60.gwf to H2-SIM-700054900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700054900-60 100% 3709KB 148.4KB/s   00:25    
Fetching /home/tania/Popcorn/data2/H2-SIM-700054960-60.gwf to H2-SIM-700054960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700054960-60 100% 3708KB 927.1KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700055020-60.gwf to H2-SIM-700055020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700055020-60 100% 3705KB 617.6KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700055080-60.gwf to H2-SIM-700055080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700055080-60 100% 3707KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700055140-60.gwf to H2-SIM-700055140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700055140-60 100% 3712KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700055200-60.gwf to H2-SIM-700055200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700055200-60 100% 3709KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700055260-60.gwf to H2-SIM-700055260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700055260-60 100% 3710KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700055320-60.gwf to H2-SIM-700055320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700055320-60 100% 3718KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700055380-60.gwf to H2-SIM-700055380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700055380-60 100% 3712KB 927.9KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700055440-60.gwf to H2-SIM-700055440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700055440-60 100% 3718KB 743.6KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700055500-60.gwf to H2-SIM-700055500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700055500-60 100% 3718KB 743.7KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700055560-60.gwf to H2-SIM-700055560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700055560-60 100% 3704KB 529.2KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700055620-60.gwf to H2-SIM-700055620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700055620-60 100% 3716KB 123.9KB/s   00:30    
Fetching /home/tania/Popcorn/data2/H2-SIM-700055680-60.gwf to H2-SIM-700055680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700055680-60 100% 3707KB 741.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700055740-60.gwf to H2-SIM-700055740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700055740-60 100% 3707KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700055800-60.gwf to H2-SIM-700055800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700055800-60 100% 3692KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700055860-60.gwf to H2-SIM-700055860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700055860-60 100% 3707KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700055920-60.gwf to H2-SIM-700055920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700055920-60 100% 3711KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700055980-60.gwf to H2-SIM-700055980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700055980-60 100% 3704KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700056040-60.gwf to H2-SIM-700056040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700056040-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700056100-60.gwf to H2-SIM-700056100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700056100-60 100% 3706KB 617.6KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700056160-60.gwf to H2-SIM-700056160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700056160-60 100% 3714KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700056220-60.gwf to H2-SIM-700056220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700056220-60 100% 3710KB 927.5KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700056280-60.gwf to H2-SIM-700056280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700056280-60 100% 3700KB 370.0KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700056340-60.gwf to H2-SIM-700056340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700056340-60 100% 3713KB 137.5KB/s   00:27    
Fetching /home/tania/Popcorn/data2/H2-SIM-700056400-60.gwf to H2-SIM-700056400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700056400-60 100% 3723KB 744.6KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700056460-60.gwf to H2-SIM-700056460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700056460-60 100% 3704KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700056520-60.gwf to H2-SIM-700056520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700056520-60 100% 3711KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700056580-60.gwf to H2-SIM-700056580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700056580-60 100% 3715KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700056640-60.gwf to H2-SIM-700056640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700056640-60 100% 3711KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700056700-60.gwf to H2-SIM-700056700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700056700-60 100% 3711KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700056760-60.gwf to H2-SIM-700056760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700056760-60 100% 3712KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700056820-60.gwf to H2-SIM-700056820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700056820-60 100% 3702KB 925.5KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700056880-60.gwf to H2-SIM-700056880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700056880-60 100% 3723KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700056940-60.gwf to H2-SIM-700056940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700056940-60 100% 3712KB 742.4KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700057000-60.gwf to H2-SIM-700057000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700057000-60 100% 3717KB 413.0KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700057060-60.gwf to H2-SIM-700057060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700057060-60 100% 3698KB 119.3KB/s   00:31    
Fetching /home/tania/Popcorn/data2/H2-SIM-700057120-60.gwf to H2-SIM-700057120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700057120-60 100% 3701KB 616.8KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700057180-60.gwf to H2-SIM-700057180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700057180-60 100% 3707KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700057240-60.gwf to H2-SIM-700057240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700057240-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700057300-60.gwf to H2-SIM-700057300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700057300-60 100% 3708KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700057360-60.gwf to H2-SIM-700057360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700057360-60 100% 3701KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700057420-60.gwf to H2-SIM-700057420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700057420-60 100% 3704KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700057480-60.gwf to H2-SIM-700057480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700057480-60 100% 3709KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700057540-60.gwf to H2-SIM-700057540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700057540-60 100% 3710KB 741.9KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700057600-60.gwf to H2-SIM-700057600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700057600-60 100% 3701KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700057660-60.gwf to H2-SIM-700057660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700057660-60 100% 3711KB 927.8KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700057720-60.gwf to H2-SIM-700057720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700057720-60 100% 3696KB 264.0KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700057780-60.gwf to H2-SIM-700057780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700057780-60 100% 3705KB 119.5KB/s   00:31    
Fetching /home/tania/Popcorn/data2/H2-SIM-700057840-60.gwf to H2-SIM-700057840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700057840-60 100% 3702KB 740.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700057900-60.gwf to H2-SIM-700057900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700057900-60 100% 3698KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700057960-60.gwf to H2-SIM-700057960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700057960-60 100% 3711KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700058020-60.gwf to H2-SIM-700058020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700058020-60 100% 3714KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700058080-60.gwf to H2-SIM-700058080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700058080-60 100% 3712KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700058140-60.gwf to H2-SIM-700058140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700058140-60 100% 3702KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700058200-60.gwf to H2-SIM-700058200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700058200-60 100% 3701KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700058260-60.gwf to H2-SIM-700058260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700058260-60 100% 3703KB 925.9KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700058320-60.gwf to H2-SIM-700058320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700058320-60 100% 3711KB 927.7KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700058380-60.gwf to H2-SIM-700058380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700058380-60 100% 3702KB 740.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700058440-60.gwf to H2-SIM-700058440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700058440-60 100% 3716KB 371.6KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700058500-60.gwf to H2-SIM-700058500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700058500-60 100% 3712KB 128.0KB/s   00:29    
Fetching /home/tania/Popcorn/data2/H2-SIM-700058560-60.gwf to H2-SIM-700058560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700058560-60 100% 3718KB 929.6KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700058620-60.gwf to H2-SIM-700058620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700058620-60 100% 3705KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700058680-60.gwf to H2-SIM-700058680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700058680-60 100% 3706KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700058740-60.gwf to H2-SIM-700058740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700058740-60 100% 3711KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700058800-60.gwf to H2-SIM-700058800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700058800-60 100% 3716KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700058860-60.gwf to H2-SIM-700058860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700058860-60 100% 3716KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700058920-60.gwf to H2-SIM-700058920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700058920-60 100% 3708KB 926.9KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700058980-60.gwf to H2-SIM-700058980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700058980-60 100% 3721KB 930.2KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700059040-60.gwf to H2-SIM-700059040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700059040-60 100% 3708KB 927.0KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700059100-60.gwf to H2-SIM-700059100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700059100-60 100% 3710KB 618.3KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700059160-60.gwf to H2-SIM-700059160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700059160-60 100% 3705KB 336.8KB/s   00:11    
Fetching /home/tania/Popcorn/data2/H2-SIM-700059220-60.gwf to H2-SIM-700059220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700059220-60 100% 3710KB 119.7KB/s   00:31    
Fetching /home/tania/Popcorn/data2/H2-SIM-700059280-60.gwf to H2-SIM-700059280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700059280-60 100% 3712KB 742.4KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700059340-60.gwf to H2-SIM-700059340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700059340-60 100% 3707KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700059400-60.gwf to H2-SIM-700059400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700059400-60 100% 3707KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700059460-60.gwf to H2-SIM-700059460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700059460-60 100% 3706KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700059520-60.gwf to H2-SIM-700059520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700059520-60 100% 3701KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700059580-60.gwf to H2-SIM-700059580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700059580-60 100% 3714KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700059640-60.gwf to H2-SIM-700059640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700059640-60 100% 3716KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700059700-60.gwf to H2-SIM-700059700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700059700-60 100% 3709KB 741.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700059760-60.gwf to H2-SIM-700059760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700059760-60 100% 3708KB 927.0KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700059820-60.gwf to H2-SIM-700059820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700059820-60 100% 3708KB 618.0KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700059880-60.gwf to H2-SIM-700059880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700059880-60 100% 3698KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700059940-60.gwf to H2-SIM-700059940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700059940-60 100% 3713KB 464.1KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700060000-60.gwf to H2-SIM-700060000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700060000-60 100% 3713KB 116.0KB/s   00:32    
Fetching /home/tania/Popcorn/data2/H2-SIM-700060060-60.gwf to H2-SIM-700060060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700060060-60 100% 3715KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700060120-60.gwf to H2-SIM-700060120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700060120-60 100% 3713KB 928.3KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700060180-60.gwf to H2-SIM-700060180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700060180-60 100% 3706KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700060240-60.gwf to H2-SIM-700060240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700060240-60 100% 3705KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700060300-60.gwf to H2-SIM-700060300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700060300-60 100% 3704KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700060360-60.gwf to H2-SIM-700060360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700060360-60 100% 3708KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700060420-60.gwf to H2-SIM-700060420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700060420-60 100% 3715KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700060480-60.gwf to H2-SIM-700060480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700060480-60 100% 3710KB 927.5KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700060540-60.gwf to H2-SIM-700060540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700060540-60 100% 3708KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700060600-60.gwf to H2-SIM-700060600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700060600-60 100% 3711KB 618.5KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700060660-60.gwf to H2-SIM-700060660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700060660-60 100% 3710KB 741.9KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700060720-60.gwf to H2-SIM-700060720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700060720-60 100% 3703KB 925.8KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700060780-60.gwf to H2-SIM-700060780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700060780-60 100% 3714KB 154.8KB/s   00:24    
Fetching /home/tania/Popcorn/data2/H2-SIM-700060840-60.gwf to H2-SIM-700060840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700060840-60 100% 3713KB 742.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700060900-60.gwf to H2-SIM-700060900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700060900-60 100% 3704KB 926.1KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700060960-60.gwf to H2-SIM-700060960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700060960-60 100% 3698KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700061020-60.gwf to H2-SIM-700061020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700061020-60 100% 3709KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700061080-60.gwf to H2-SIM-700061080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700061080-60 100% 3717KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700061140-60.gwf to H2-SIM-700061140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700061140-60 100% 3721KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700061200-60.gwf to H2-SIM-700061200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700061200-60 100% 3707KB 741.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700061260-60.gwf to H2-SIM-700061260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700061260-60 100% 3705KB 926.3KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700061320-60.gwf to H2-SIM-700061320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700061320-60 100% 3707KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700061380-60.gwf to H2-SIM-700061380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700061380-60 100% 3709KB 741.7KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700061440-60.gwf to H2-SIM-700061440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700061440-60 100% 3690KB 263.6KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700061500-60.gwf to H2-SIM-700061500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700061500-60 100% 3718KB 148.7KB/s   00:25    
Fetching /home/tania/Popcorn/data2/H2-SIM-700061560-60.gwf to H2-SIM-700061560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700061560-60 100% 3706KB 617.7KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700061620-60.gwf to H2-SIM-700061620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700061620-60 100% 3706KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700061680-60.gwf to H2-SIM-700061680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700061680-60 100% 3715KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700061740-60.gwf to H2-SIM-700061740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700061740-60 100% 3710KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700061800-60.gwf to H2-SIM-700061800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700061800-60 100% 3715KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700061860-60.gwf to H2-SIM-700061860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700061860-60 100% 3710KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700061920-60.gwf to H2-SIM-700061920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700061920-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700061980-60.gwf to H2-SIM-700061980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700061980-60 100% 3703KB 740.6KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700062040-60.gwf to H2-SIM-700062040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700062040-60 100% 3709KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700062100-60.gwf to H2-SIM-700062100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700062100-60 100% 3714KB 742.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700062160-60.gwf to H2-SIM-700062160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700062160-60 100% 3706KB 926.4KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700062220-60.gwf to H2-SIM-700062220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700062220-60 100% 3711KB 530.2KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700062280-60.gwf to H2-SIM-700062280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700062280-60 100% 3718KB 120.0KB/s   00:31    
Fetching /home/tania/Popcorn/data2/H2-SIM-700062340-60.gwf to H2-SIM-700062340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700062340-60 100% 3715KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700062400-60.gwf to H2-SIM-700062400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700062400-60 100% 3692KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700062460-60.gwf to H2-SIM-700062460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700062460-60 100% 3712KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700062520-60.gwf to H2-SIM-700062520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700062520-60 100% 3715KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700062580-60.gwf to H2-SIM-700062580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700062580-60 100% 3712KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700062640-60.gwf to H2-SIM-700062640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700062640-60 100% 3716KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700062700-60.gwf to H2-SIM-700062700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700062700-60 100% 3703KB 370.3KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700062760-60.gwf to H2-SIM-700062760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700062760-60 100% 3710KB 142.7KB/s   00:26    
Fetching /home/tania/Popcorn/data2/H2-SIM-700062820-60.gwf to H2-SIM-700062820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700062820-60 100% 3716KB 743.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700062880-60.gwf to H2-SIM-700062880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700062880-60 100% 3713KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700062940-60.gwf to H2-SIM-700062940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700062940-60 100% 3719KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700063000-60.gwf to H2-SIM-700063000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700063000-60 100% 3720KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700063060-60.gwf to H2-SIM-700063060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700063060-60 100% 3709KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700063120-60.gwf to H2-SIM-700063120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700063120-60 100% 3713KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700063180-60.gwf to H2-SIM-700063180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700063180-60 100% 3714KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700063240-60.gwf to H2-SIM-700063240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700063240-60 100% 3712KB 928.0KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700063300-60.gwf to H2-SIM-700063300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700063300-60 100% 3699KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700063360-60.gwf to H2-SIM-700063360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700063360-60 100% 3704KB 740.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700063420-60.gwf to H2-SIM-700063420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700063420-60 100% 3704KB 370.4KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700063480-60.gwf to H2-SIM-700063480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700063480-60 100% 3697KB 127.5KB/s   00:29    
Fetching /home/tania/Popcorn/data2/H2-SIM-700063540-60.gwf to H2-SIM-700063540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700063540-60 100% 3709KB 618.2KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700063600-60.gwf to H2-SIM-700063600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700063600-60 100% 3712KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700063660-60.gwf to H2-SIM-700063660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700063660-60 100% 3713KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700063720-60.gwf to H2-SIM-700063720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700063720-60 100% 3693KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700063780-60.gwf to H2-SIM-700063780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700063780-60 100% 3713KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700063840-60.gwf to H2-SIM-700063840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700063840-60 100% 3721KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700063900-60.gwf to H2-SIM-700063900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700063900-60 100% 3712KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700063960-60.gwf to H2-SIM-700063960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700063960-60 100% 3706KB 741.2KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700064020-60.gwf to H2-SIM-700064020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700064020-60 100% 3703KB 925.7KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700064080-60.gwf to H2-SIM-700064080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700064080-60 100% 3716KB 743.1KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700064140-60.gwf to H2-SIM-700064140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700064140-60 100% 3714KB 309.5KB/s   00:12    
Fetching /home/tania/Popcorn/data2/H2-SIM-700064200-60.gwf to H2-SIM-700064200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700064200-60 100% 3703KB 148.1KB/s   00:25    
Fetching /home/tania/Popcorn/data2/H2-SIM-700064260-60.gwf to H2-SIM-700064260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700064260-60 100% 3694KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700064320-60.gwf to H2-SIM-700064320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700064320-60 100% 3719KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700064380-60.gwf to H2-SIM-700064380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700064380-60 100% 3706KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700064440-60.gwf to H2-SIM-700064440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700064440-60 100% 3700KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700064500-60.gwf to H2-SIM-700064500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700064500-60 100% 3702KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700064560-60.gwf to H2-SIM-700064560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700064560-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700064620-60.gwf to H2-SIM-700064620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700064620-60 100% 3719KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700064680-60.gwf to H2-SIM-700064680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700064680-60 100% 3710KB 742.1KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700064740-60.gwf to H2-SIM-700064740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700064740-60 100% 3714KB 928.6KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700064800-60.gwf to H2-SIM-700064800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700064800-60 100% 3711KB 927.7KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700064860-60.gwf to H2-SIM-700064860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700064860-60 100% 3715KB 371.5KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700064920-60.gwf to H2-SIM-700064920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700064920-60 100% 3705KB 119.5KB/s   00:31    
Fetching /home/tania/Popcorn/data2/H2-SIM-700064980-60.gwf to H2-SIM-700064980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700064980-60 100% 3711KB 742.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700065040-60.gwf to H2-SIM-700065040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700065040-60 100% 3697KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700065100-60.gwf to H2-SIM-700065100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700065100-60 100% 3697KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700065160-60.gwf to H2-SIM-700065160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700065160-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700065220-60.gwf to H2-SIM-700065220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700065220-60 100% 3710KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700065280-60.gwf to H2-SIM-700065280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700065280-60 100% 3710KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700065340-60.gwf to H2-SIM-700065340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700065340-60 100% 3719KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700065400-60.gwf to H2-SIM-700065400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700065400-60 100% 3700KB 411.1KB/s   00:09    
Fetching /home/tania/Popcorn/data2/H2-SIM-700065460-60.gwf to H2-SIM-700065460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700065460-60 100% 3709KB 176.6KB/s   00:21    
Fetching /home/tania/Popcorn/data2/H2-SIM-700065520-60.gwf to H2-SIM-700065520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700065520-60 100% 3713KB 285.6KB/s   00:13    
Fetching /home/tania/Popcorn/data2/H2-SIM-700065580-60.gwf to H2-SIM-700065580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700065580-60 100% 3708KB 741.6KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700065640-60.gwf to H2-SIM-700065640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700065640-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700065700-60.gwf to H2-SIM-700065700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700065700-60 100% 3710KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700065760-60.gwf to H2-SIM-700065760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700065760-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700065820-60.gwf to H2-SIM-700065820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700065820-60 100% 3711KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700065880-60.gwf to H2-SIM-700065880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700065880-60 100% 3702KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700065940-60.gwf to H2-SIM-700065940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700065940-60 100% 3707KB 741.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700066000-60.gwf to H2-SIM-700066000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700066000-60 100% 3710KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700066060-60.gwf to H2-SIM-700066060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700066060-60 100% 3713KB 337.6KB/s   00:11    
Fetching /home/tania/Popcorn/data2/H2-SIM-700066120-60.gwf to H2-SIM-700066120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700066120-60 100% 3698KB 108.8KB/s   00:34    
Fetching /home/tania/Popcorn/data2/H2-SIM-700066180-60.gwf to H2-SIM-700066180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700066180-60 100% 3711KB 927.7KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700066240-60.gwf to H2-SIM-700066240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700066240-60 100% 3711KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700066300-60.gwf to H2-SIM-700066300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700066300-60 100% 3712KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700066360-60.gwf to H2-SIM-700066360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700066360-60 100% 3719KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700066420-60.gwf to H2-SIM-700066420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700066420-60 100% 3711KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700066480-60.gwf to H2-SIM-700066480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700066480-60 100% 3700KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700066540-60.gwf to H2-SIM-700066540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700066540-60 100% 3704KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700066600-60.gwf to H2-SIM-700066600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700066600-60 100% 3700KB 925.1KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700066660-60.gwf to H2-SIM-700066660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700066660-60 100% 3714KB 928.5KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700066720-60.gwf to H2-SIM-700066720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700066720-60 100% 3710KB 742.0KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700066780-60.gwf to H2-SIM-700066780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700066780-60 100% 3717KB 195.6KB/s   00:19    
Fetching /home/tania/Popcorn/data2/H2-SIM-700066840-60.gwf to H2-SIM-700066840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700066840-60 100% 3715KB 154.8KB/s   00:24    
Fetching /home/tania/Popcorn/data2/H2-SIM-700066900-60.gwf to H2-SIM-700066900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700066900-60 100% 3712KB 742.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700066960-60.gwf to H2-SIM-700066960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700066960-60 100% 3711KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700067020-60.gwf to H2-SIM-700067020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700067020-60 100% 3704KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700067080-60.gwf to H2-SIM-700067080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700067080-60 100% 3702KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700067140-60.gwf to H2-SIM-700067140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700067140-60 100% 3704KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700067200-60.gwf to H2-SIM-700067200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700067200-60 100% 3720KB 744.0KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700067260-60.gwf to H2-SIM-700067260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700067260-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700067320-60.gwf to H2-SIM-700067320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700067320-60 100% 3697KB 739.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700067380-60.gwf to H2-SIM-700067380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700067380-60 100% 3709KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700067440-60.gwf to H2-SIM-700067440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700067440-60 100% 3714KB 371.4KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700067500-60.gwf to H2-SIM-700067500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700067500-60 100% 3727KB 133.1KB/s   00:28    
Fetching /home/tania/Popcorn/data2/H2-SIM-700067560-60.gwf to H2-SIM-700067560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700067560-60 100% 3707KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700067620-60.gwf to H2-SIM-700067620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700067620-60 100% 3714KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700067680-60.gwf to H2-SIM-700067680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700067680-60 100% 3699KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700067740-60.gwf to H2-SIM-700067740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700067740-60 100% 3706KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700067800-60.gwf to H2-SIM-700067800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700067800-60 100% 3710KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700067860-60.gwf to H2-SIM-700067860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700067860-60 100% 3714KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700067920-60.gwf to H2-SIM-700067920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700067920-60 100% 3715KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700067980-60.gwf to H2-SIM-700067980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700067980-60 100% 3718KB 743.6KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700068040-60.gwf to H2-SIM-700068040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700068040-60 100% 3718KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700068100-60.gwf to H2-SIM-700068100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700068100-60 100% 3709KB 741.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700068160-60.gwf to H2-SIM-700068160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700068160-60 100% 3710KB 741.9KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700068220-60.gwf to H2-SIM-700068220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700068220-60 100% 3708KB 337.1KB/s   00:11    
Fetching /home/tania/Popcorn/data2/H2-SIM-700068280-60.gwf to H2-SIM-700068280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700068280-60 100% 3723KB 124.1KB/s   00:30    
Fetching /home/tania/Popcorn/data2/H2-SIM-700068340-60.gwf to H2-SIM-700068340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700068340-60 100% 3699KB 924.8KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700068400-60.gwf to H2-SIM-700068400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700068400-60 100% 3715KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700068460-60.gwf to H2-SIM-700068460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700068460-60 100% 3707KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700068520-60.gwf to H2-SIM-700068520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700068520-60 100% 3712KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700068580-60.gwf to H2-SIM-700068580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700068580-60 100% 3715KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700068640-60.gwf to H2-SIM-700068640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700068640-60 100% 3717KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700068700-60.gwf to H2-SIM-700068700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700068700-60 100% 3696KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700068760-60.gwf to H2-SIM-700068760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700068760-60 100% 3710KB 618.4KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700068820-60.gwf to H2-SIM-700068820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700068820-60 100% 3707KB 926.8KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700068880-60.gwf to H2-SIM-700068880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700068880-60 100% 3711KB 742.2KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700068940-60.gwf to H2-SIM-700068940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700068940-60 100% 3713KB 928.2KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700069000-60.gwf to H2-SIM-700069000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700069000-60 100% 3707KB 529.5KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700069060-60.gwf to H2-SIM-700069060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700069060-60 100% 3695KB 105.6KB/s   00:35    
Fetching /home/tania/Popcorn/data2/H2-SIM-700069120-60.gwf to H2-SIM-700069120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700069120-60 100% 3709KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700069180-60.gwf to H2-SIM-700069180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700069180-60 100% 3706KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700069240-60.gwf to H2-SIM-700069240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700069240-60 100% 3712KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700069300-60.gwf to H2-SIM-700069300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700069300-60 100% 3706KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700069360-60.gwf to H2-SIM-700069360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700069360-60 100% 3709KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700069420-60.gwf to H2-SIM-700069420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700069420-60 100% 3724KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700069480-60.gwf to H2-SIM-700069480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700069480-60 100% 3718KB 619.7KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700069540-60.gwf to H2-SIM-700069540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700069540-60 100% 3715KB 185.8KB/s   00:20    
Fetching /home/tania/Popcorn/data2/H2-SIM-700069600-60.gwf to H2-SIM-700069600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700069600-60 100% 3700KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700069660-60.gwf to H2-SIM-700069660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700069660-60 100% 3719KB 929.7KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700069720-60.gwf to H2-SIM-700069720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700069720-60 100% 3711KB 742.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700069780-60.gwf to H2-SIM-700069780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700069780-60 100% 3711KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700069840-60.gwf to H2-SIM-700069840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700069840-60 100% 3705KB 617.5KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700069900-60.gwf to H2-SIM-700069900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700069900-60 100% 3695KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700069960-60.gwf to H2-SIM-700069960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700069960-60 100% 3702KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700070020-60.gwf to H2-SIM-700070020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700070020-60 100% 3705KB 741.0KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700070080-60.gwf to H2-SIM-700070080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700070080-60 100% 3713KB 618.9KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700070140-60.gwf to H2-SIM-700070140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700070140-60 100% 3721KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700070200-60.gwf to H2-SIM-700070200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700070200-60 100% 3708KB 927.1KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700070260-60.gwf to H2-SIM-700070260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700070260-60 100% 3701KB 370.1KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700070320-60.gwf to H2-SIM-700070320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700070320-60 100% 3704KB 127.7KB/s   00:29    
Fetching /home/tania/Popcorn/data2/H2-SIM-700070380-60.gwf to H2-SIM-700070380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700070380-60 100% 3715KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700070440-60.gwf to H2-SIM-700070440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700070440-60 100% 3717KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700070500-60.gwf to H2-SIM-700070500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700070500-60 100% 3720KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700070560-60.gwf to H2-SIM-700070560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700070560-60 100% 3709KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700070620-60.gwf to H2-SIM-700070620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700070620-60 100% 3715KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700070680-60.gwf to H2-SIM-700070680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700070680-60 100% 3703KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700070740-60.gwf to H2-SIM-700070740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700070740-60 100% 3720KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700070800-60.gwf to H2-SIM-700070800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700070800-60 100% 3706KB 926.6KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700070860-60.gwf to H2-SIM-700070860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700070860-60 100% 3706KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700070920-60.gwf to H2-SIM-700070920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700070920-60 100% 3720KB 743.9KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700070980-60.gwf to H2-SIM-700070980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700070980-60 100% 3716KB 619.4KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700071040-60.gwf to H2-SIM-700071040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700071040-60 100% 3698KB 739.7KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700071100-60.gwf to H2-SIM-700071100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700071100-60 100% 3711KB 247.4KB/s   00:15    
Fetching /home/tania/Popcorn/data2/H2-SIM-700071160-60.gwf to H2-SIM-700071160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700071160-60 100% 3699KB 739.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700071220-60.gwf to H2-SIM-700071220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700071220-60 100% 3710KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700071280-60.gwf to H2-SIM-700071280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700071280-60 100% 3704KB 740.9KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700071340-60.gwf to H2-SIM-700071340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700071340-60 100% 3712KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700071400-60.gwf to H2-SIM-700071400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700071400-60 100% 3718KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700071460-60.gwf to H2-SIM-700071460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700071460-60 100% 3707KB 463.4KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700071520-60.gwf to H2-SIM-700071520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700071520-60 100% 3710KB 309.2KB/s   00:12    
Fetching /home/tania/Popcorn/data2/H2-SIM-700071580-60.gwf to H2-SIM-700071580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700071580-60 100% 3711KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700071640-60.gwf to H2-SIM-700071640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700071640-60 100% 3715KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700071700-60.gwf to H2-SIM-700071700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700071700-60 100% 3707KB 617.8KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700071760-60.gwf to H2-SIM-700071760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700071760-60 100% 3708KB 741.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700071820-60.gwf to H2-SIM-700071820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700071820-60 100% 3708KB 154.5KB/s   00:24    
Fetching /home/tania/Popcorn/data2/H2-SIM-700071880-60.gwf to H2-SIM-700071880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700071880-60 100% 3705KB 741.1KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700071940-60.gwf to H2-SIM-700071940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700071940-60 100% 3709KB 741.9KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700072000-60.gwf to H2-SIM-700072000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700072000-60 100% 3711KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700072060-60.gwf to H2-SIM-700072060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700072060-60 100% 3716KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700072120-60.gwf to H2-SIM-700072120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700072120-60 100% 3714KB 928.5KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700072180-60.gwf to H2-SIM-700072180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700072180-60 100% 3717KB 929.3KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700072240-60.gwf to H2-SIM-700072240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700072240-60 100% 3717KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700072300-60.gwf to H2-SIM-700072300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700072300-60 100% 3714KB 928.4KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700072360-60.gwf to H2-SIM-700072360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700072360-60 100% 3723KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700072420-60.gwf to H2-SIM-700072420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700072420-60 100% 3713KB 618.8KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700072480-60.gwf to H2-SIM-700072480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700072480-60 100% 3709KB 154.6KB/s   00:24    
Fetching /home/tania/Popcorn/data2/H2-SIM-700072540-60.gwf to H2-SIM-700072540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700072540-60 100% 3715KB 742.9KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700072600-60.gwf to H2-SIM-700072600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700072600-60 100% 3709KB 927.3KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700072660-60.gwf to H2-SIM-700072660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700072660-60 100% 3713KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700072720-60.gwf to H2-SIM-700072720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700072720-60 100% 3710KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700072780-60.gwf to H2-SIM-700072780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700072780-60 100% 3713KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700072840-60.gwf to H2-SIM-700072840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700072840-60 100% 3700KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700072900-60.gwf to H2-SIM-700072900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700072900-60 100% 3724KB 744.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700072960-60.gwf to H2-SIM-700072960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700072960-60 100% 3697KB 739.4KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700073020-60.gwf to H2-SIM-700073020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700073020-60 100% 3701KB 925.1KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700073080-60.gwf to H2-SIM-700073080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700073080-60 100% 3716KB 928.9KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700073140-60.gwf to H2-SIM-700073140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700073140-60 100% 3719KB 464.9KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700073200-60.gwf to H2-SIM-700073200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700073200-60 100% 3704KB 137.2KB/s   00:27    
Fetching /home/tania/Popcorn/data2/H2-SIM-700073260-60.gwf to H2-SIM-700073260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700073260-60 100% 3720KB 744.0KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700073320-60.gwf to H2-SIM-700073320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700073320-60 100% 3704KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700073380-60.gwf to H2-SIM-700073380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700073380-60 100% 3712KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700073440-60.gwf to H2-SIM-700073440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700073440-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700073500-60.gwf to H2-SIM-700073500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700073500-60 100% 3707KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700073560-60.gwf to H2-SIM-700073560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700073560-60 100% 3720KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700073620-60.gwf to H2-SIM-700073620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700073620-60 100% 3717KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700073680-60.gwf to H2-SIM-700073680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700073680-60 100% 3704KB 926.0KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700073740-60.gwf to H2-SIM-700073740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700073740-60 100% 3721KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700073800-60.gwf to H2-SIM-700073800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700073800-60 100% 3722KB 744.4KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700073860-60.gwf to H2-SIM-700073860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700073860-60 100% 3717KB 743.4KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700073920-60.gwf to H2-SIM-700073920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700073920-60 100% 3705KB 529.2KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700073980-60.gwf to H2-SIM-700073980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700073980-60 100% 3703KB 119.5KB/s   00:31    
Fetching /home/tania/Popcorn/data2/H2-SIM-700074040-60.gwf to H2-SIM-700074040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700074040-60 100% 3710KB 618.3KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700074100-60.gwf to H2-SIM-700074100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700074100-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700074160-60.gwf to H2-SIM-700074160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700074160-60 100% 3703KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700074220-60.gwf to H2-SIM-700074220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700074220-60 100% 3700KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700074280-60.gwf to H2-SIM-700074280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700074280-60 100% 3720KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700074340-60.gwf to H2-SIM-700074340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700074340-60 100% 3709KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700074400-60.gwf to H2-SIM-700074400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700074400-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700074460-60.gwf to H2-SIM-700074460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700074460-60 100% 3701KB 616.8KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700074520-60.gwf to H2-SIM-700074520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700074520-60 100% 3705KB 926.2KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700074580-60.gwf to H2-SIM-700074580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700074580-60 100% 3713KB 928.3KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700074640-60.gwf to H2-SIM-700074640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700074640-60 100% 3718KB 286.0KB/s   00:13    
Fetching /home/tania/Popcorn/data2/H2-SIM-700074700-60.gwf to H2-SIM-700074700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700074700-60 100% 3714KB 116.1KB/s   00:32    
Fetching /home/tania/Popcorn/data2/H2-SIM-700074760-60.gwf to H2-SIM-700074760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700074760-60 100% 3696KB 739.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700074820-60.gwf to H2-SIM-700074820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700074820-60 100% 3723KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700074880-60.gwf to H2-SIM-700074880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700074880-60 100% 3710KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700074940-60.gwf to H2-SIM-700074940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700074940-60 100% 3717KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700075000-60.gwf to H2-SIM-700075000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700075000-60 100% 3711KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700075060-60.gwf to H2-SIM-700075060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700075060-60 100% 3710KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700075120-60.gwf to H2-SIM-700075120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700075120-60 100% 3709KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700075180-60.gwf to H2-SIM-700075180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700075180-60 100% 3710KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700075240-60.gwf to H2-SIM-700075240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700075240-60 100% 3713KB 928.3KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700075300-60.gwf to H2-SIM-700075300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700075300-60 100% 3702KB 925.4KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700075360-60.gwf to H2-SIM-700075360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700075360-60 100% 3711KB 247.4KB/s   00:15    
Fetching /home/tania/Popcorn/data2/H2-SIM-700075420-60.gwf to H2-SIM-700075420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700075420-60 100% 3693KB 147.7KB/s   00:25    
Fetching /home/tania/Popcorn/data2/H2-SIM-700075480-60.gwf to H2-SIM-700075480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700075480-60 100% 3704KB 617.4KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700075540-60.gwf to H2-SIM-700075540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700075540-60 100% 3709KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700075600-60.gwf to H2-SIM-700075600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700075600-60 100% 3719KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700075660-60.gwf to H2-SIM-700075660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700075660-60 100% 3720KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700075720-60.gwf to H2-SIM-700075720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700075720-60 100% 3713KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700075780-60.gwf to H2-SIM-700075780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700075780-60 100% 3716KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700075840-60.gwf to H2-SIM-700075840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700075840-60 100% 3710KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700075900-60.gwf to H2-SIM-700075900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700075900-60 100% 3715KB 619.2KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700075960-60.gwf to H2-SIM-700075960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700075960-60 100% 3718KB 929.6KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700076020-60.gwf to H2-SIM-700076020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700076020-60 100% 3708KB 617.9KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700076080-60.gwf to H2-SIM-700076080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700076080-60 100% 3717KB 743.4KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700076140-60.gwf to H2-SIM-700076140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700076140-60 100% 3705KB 529.3KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700076200-60.gwf to H2-SIM-700076200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700076200-60 100% 3720KB 124.0KB/s   00:30    
Fetching /home/tania/Popcorn/data2/H2-SIM-700076260-60.gwf to H2-SIM-700076260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700076260-60 100% 3707KB 741.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700076320-60.gwf to H2-SIM-700076320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700076320-60 100% 3717KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700076380-60.gwf to H2-SIM-700076380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700076380-60 100% 3706KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700076440-60.gwf to H2-SIM-700076440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700076440-60 100% 3717KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700076500-60.gwf to H2-SIM-700076500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700076500-60 100% 3714KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700076560-60.gwf to H2-SIM-700076560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700076560-60 100% 3703KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700076620-60.gwf to H2-SIM-700076620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700076620-60 100% 3706KB 529.4KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700076680-60.gwf to H2-SIM-700076680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700076680-60 100% 3721KB 161.8KB/s   00:23    
Fetching /home/tania/Popcorn/data2/H2-SIM-700076740-60.gwf to H2-SIM-700076740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700076740-60 100% 3702KB 740.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700076800-60.gwf to H2-SIM-700076800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700076800-60 100% 3720KB 620.0KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700076860-60.gwf to H2-SIM-700076860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700076860-60 100% 3705KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700076920-60.gwf to H2-SIM-700076920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700076920-60 100% 3714KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700076980-60.gwf to H2-SIM-700076980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700076980-60 100% 3714KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700077040-60.gwf to H2-SIM-700077040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700077040-60 100% 3708KB 741.6KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700077100-60.gwf to H2-SIM-700077100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700077100-60 100% 3717KB 743.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700077160-60.gwf to H2-SIM-700077160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700077160-60 100% 3715KB 928.9KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700077220-60.gwf to H2-SIM-700077220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700077220-60 100% 3713KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700077280-60.gwf to H2-SIM-700077280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700077280-60 100% 3706KB 926.6KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700077340-60.gwf to H2-SIM-700077340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700077340-60 100% 3698KB 616.3KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700077400-60.gwf to H2-SIM-700077400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700077400-60 100% 3708KB 741.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700077460-60.gwf to H2-SIM-700077460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700077460-60 100% 3709KB 927.2KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700077520-60.gwf to H2-SIM-700077520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700077520-60 100% 3699KB 616.6KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700077580-60.gwf to H2-SIM-700077580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700077580-60 100% 3721KB 744.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700077640-60.gwf to H2-SIM-700077640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700077640-60 100% 3708KB 206.0KB/s   00:18    
Fetching /home/tania/Popcorn/data2/H2-SIM-700077700-60.gwf to H2-SIM-700077700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700077700-60 100% 3698KB 924.5KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700077760-60.gwf to H2-SIM-700077760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700077760-60 100% 3699KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700077820-60.gwf to H2-SIM-700077820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700077820-60 100% 3706KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700077880-60.gwf to H2-SIM-700077880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700077880-60 100% 3698KB 336.2KB/s   00:11    
Fetching /home/tania/Popcorn/data2/H2-SIM-700077940-60.gwf to H2-SIM-700077940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700077940-60 100% 3709KB 119.6KB/s   00:31    
Fetching /home/tania/Popcorn/data2/H2-SIM-700078000-60.gwf to H2-SIM-700078000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700078000-60 100% 3708KB 927.0KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700078060-60.gwf to H2-SIM-700078060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700078060-60 100% 3703KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700078120-60.gwf to H2-SIM-700078120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700078120-60 100% 3711KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700078180-60.gwf to H2-SIM-700078180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700078180-60 100% 3709KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700078240-60.gwf to H2-SIM-700078240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700078240-60 100% 3702KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700078300-60.gwf to H2-SIM-700078300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700078300-60 100% 3703KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700078360-60.gwf to H2-SIM-700078360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700078360-60 100% 3704KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700078420-60.gwf to H2-SIM-700078420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700078420-60 100% 3705KB 741.1KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700078480-60.gwf to H2-SIM-700078480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700078480-60 100% 3721KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700078540-60.gwf to H2-SIM-700078540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700078540-60 100% 3703KB 925.8KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700078600-60.gwf to H2-SIM-700078600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700078600-60 100% 3708KB 309.0KB/s   00:12    
Fetching /home/tania/Popcorn/data2/H2-SIM-700078660-60.gwf to H2-SIM-700078660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700078660-60 100% 3705KB 127.8KB/s   00:29    
Fetching /home/tania/Popcorn/data2/H2-SIM-700078720-60.gwf to H2-SIM-700078720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700078720-60 100% 3709KB 927.2KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700078780-60.gwf to H2-SIM-700078780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700078780-60 100% 3701KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700078840-60.gwf to H2-SIM-700078840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700078840-60 100% 3713KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700078900-60.gwf to H2-SIM-700078900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700078900-60 100% 3708KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700078960-60.gwf to H2-SIM-700078960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700078960-60 100% 3705KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700079020-60.gwf to H2-SIM-700079020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700079020-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700079080-60.gwf to H2-SIM-700079080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700079080-60 100% 3714KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700079140-60.gwf to H2-SIM-700079140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700079140-60 100% 3707KB 370.7KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700079200-60.gwf to H2-SIM-700079200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700079200-60 100% 3716KB 137.6KB/s   00:27    
Fetching /home/tania/Popcorn/data2/H2-SIM-700079260-60.gwf to H2-SIM-700079260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700079260-60 100% 3698KB 528.2KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700079320-60.gwf to H2-SIM-700079320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700079320-60 100% 3715KB 742.9KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700079380-60.gwf to H2-SIM-700079380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700079380-60 100% 3709KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700079440-60.gwf to H2-SIM-700079440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700079440-60 100% 3710KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700079500-60.gwf to H2-SIM-700079500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700079500-60 100% 3713KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700079560-60.gwf to H2-SIM-700079560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700079560-60 100% 3702KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700079620-60.gwf to H2-SIM-700079620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700079620-60 100% 3704KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700079680-60.gwf to H2-SIM-700079680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700079680-60 100% 3720KB 465.0KB/s   00:08    
Fetching /home/tania/Popcorn/data2/H2-SIM-700079740-60.gwf to H2-SIM-700079740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700079740-60 100% 3709KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700079800-60.gwf to H2-SIM-700079800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700079800-60 100% 3712KB 371.2KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700079860-60.gwf to H2-SIM-700079860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700079860-60 100% 3708KB 123.6KB/s   00:30    
Fetching /home/tania/Popcorn/data2/H2-SIM-700079920-60.gwf to H2-SIM-700079920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700079920-60 100% 3707KB 617.8KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700079980-60.gwf to H2-SIM-700079980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700079980-60 100% 3714KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700080040-60.gwf to H2-SIM-700080040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700080040-60 100% 3706KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700080100-60.gwf to H2-SIM-700080100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700080100-60 100% 3715KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700080160-60.gwf to H2-SIM-700080160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700080160-60 100% 3711KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700080220-60.gwf to H2-SIM-700080220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700080220-60 100% 3713KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700080280-60.gwf to H2-SIM-700080280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700080280-60 100% 3713KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700080340-60.gwf to H2-SIM-700080340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700080340-60 100% 3707KB 741.4KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700080400-60.gwf to H2-SIM-700080400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700080400-60 100% 3712KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700080460-60.gwf to H2-SIM-700080460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700080460-60 100% 3695KB 739.0KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700080520-60.gwf to H2-SIM-700080520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700080520-60 100% 3709KB 741.8KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700080580-60.gwf to H2-SIM-700080580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700080580-60 100% 3717KB 337.9KB/s   00:11    
Fetching /home/tania/Popcorn/data2/H2-SIM-700080640-60.gwf to H2-SIM-700080640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700080640-60 100% 3715KB 119.8KB/s   00:31    
Fetching /home/tania/Popcorn/data2/H2-SIM-700080700-60.gwf to H2-SIM-700080700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700080700-60 100% 3717KB 743.4KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700080760-60.gwf to H2-SIM-700080760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700080760-60 100% 3715KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700080820-60.gwf to H2-SIM-700080820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700080820-60 100% 3714KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700080880-60.gwf to H2-SIM-700080880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700080880-60 100% 3722KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700080940-60.gwf to H2-SIM-700080940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700080940-60 100% 3704KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700081000-60.gwf to H2-SIM-700081000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700081000-60 100% 3716KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700081060-60.gwf to H2-SIM-700081060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700081060-60 100% 3708KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700081120-60.gwf to H2-SIM-700081120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700081120-60 100% 3705KB 926.2KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700081180-60.gwf to H2-SIM-700081180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700081180-60 100% 3709KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700081240-60.gwf to H2-SIM-700081240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700081240-60 100% 3714KB 928.6KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700081300-60.gwf to H2-SIM-700081300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700081300-60 100% 3707KB 370.7KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700081360-60.gwf to H2-SIM-700081360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700081360-60 100% 3707KB 741.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700081420-60.gwf to H2-SIM-700081420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700081420-60 100% 3716KB 142.9KB/s   00:26    
Fetching /home/tania/Popcorn/data2/H2-SIM-700081480-60.gwf to H2-SIM-700081480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700081480-60 100% 3710KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700081540-60.gwf to H2-SIM-700081540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700081540-60 100% 3715KB 928.8KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700081600-60.gwf to H2-SIM-700081600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700081600-60 100% 3704KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700081660-60.gwf to H2-SIM-700081660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700081660-60 100% 3702KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700081720-60.gwf to H2-SIM-700081720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700081720-60 100% 3709KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700081780-60.gwf to H2-SIM-700081780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700081780-60 100% 3703KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700081840-60.gwf to H2-SIM-700081840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700081840-60 100% 3718KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700081900-60.gwf to H2-SIM-700081900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700081900-60 100% 3719KB 743.7KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700081960-60.gwf to H2-SIM-700081960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700081960-60 100% 3715KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700082020-60.gwf to H2-SIM-700082020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700082020-60 100% 3717KB 371.7KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700082080-60.gwf to H2-SIM-700082080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700082080-60 100% 3710KB 127.9KB/s   00:29    
Fetching /home/tania/Popcorn/data2/H2-SIM-700082140-60.gwf to H2-SIM-700082140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700082140-60 100% 3705KB 741.1KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700082200-60.gwf to H2-SIM-700082200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700082200-60 100% 3701KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700082260-60.gwf to H2-SIM-700082260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700082260-60 100% 3716KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700082320-60.gwf to H2-SIM-700082320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700082320-60 100% 3708KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700082380-60.gwf to H2-SIM-700082380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700082380-60 100% 3714KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700082440-60.gwf to H2-SIM-700082440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700082440-60 100% 3713KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700082500-60.gwf to H2-SIM-700082500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700082500-60 100% 3710KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700082560-60.gwf to H2-SIM-700082560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700082560-60 100% 3715KB 928.7KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700082620-60.gwf to H2-SIM-700082620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700082620-60 100% 3703KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700082680-60.gwf to H2-SIM-700082680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700082680-60 100% 3702KB 528.8KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700082740-60.gwf to H2-SIM-700082740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700082740-60 100% 3715KB 742.9KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700082800-60.gwf to H2-SIM-700082800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700082800-60 100% 3713KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700082860-60.gwf to H2-SIM-700082860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700082860-60 100% 3706KB 617.6KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700082920-60.gwf to H2-SIM-700082920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700082920-60 100% 3712KB 123.7KB/s   00:30    
Fetching /home/tania/Popcorn/data2/H2-SIM-700082980-60.gwf to H2-SIM-700082980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700082980-60 100% 3712KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700083040-60.gwf to H2-SIM-700083040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700083040-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700083100-60.gwf to H2-SIM-700083100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700083100-60 100% 3706KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700083160-60.gwf to H2-SIM-700083160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700083160-60 100% 3715KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700083220-60.gwf to H2-SIM-700083220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700083220-60 100% 3707KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700083280-60.gwf to H2-SIM-700083280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700083280-60 100% 3705KB 926.3KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700083340-60.gwf to H2-SIM-700083340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700083340-60 100% 3712KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700083400-60.gwf to H2-SIM-700083400-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700083400-60 100% 3708KB 741.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700083460-60.gwf to H2-SIM-700083460-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700083460-60 100% 3713KB 928.2KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700083520-60.gwf to H2-SIM-700083520-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700083520-60 100% 3706KB 741.2KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700083580-60.gwf to H2-SIM-700083580-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700083580-60 100% 3716KB 285.9KB/s   00:13    
Fetching /home/tania/Popcorn/data2/H2-SIM-700083640-60.gwf to H2-SIM-700083640-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700083640-60 100% 3707KB 123.6KB/s   00:30    
Fetching /home/tania/Popcorn/data2/H2-SIM-700083700-60.gwf to H2-SIM-700083700-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700083700-60 100% 3702KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700083760-60.gwf to H2-SIM-700083760-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700083760-60 100% 3708KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700083820-60.gwf to H2-SIM-700083820-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700083820-60 100% 3703KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700083880-60.gwf to H2-SIM-700083880-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700083880-60 100% 3710KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700083940-60.gwf to H2-SIM-700083940-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700083940-60 100% 3714KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700084000-60.gwf to H2-SIM-700084000-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700084000-60 100% 3711KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700084060-60.gwf to H2-SIM-700084060-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700084060-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700084120-60.gwf to H2-SIM-700084120-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700084120-60 100% 3702KB 740.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700084180-60.gwf to H2-SIM-700084180-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700084180-60 100% 3702KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700084240-60.gwf to H2-SIM-700084240-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700084240-60 100% 3710KB 927.5KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700084300-60.gwf to H2-SIM-700084300-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700084300-60 100% 3710KB 530.0KB/s   00:07    
Fetching /home/tania/Popcorn/data2/H2-SIM-700084360-60.gwf to H2-SIM-700084360-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700084360-60 100% 3707KB 119.6KB/s   00:31    
Fetching /home/tania/Popcorn/data2/H2-SIM-700084420-60.gwf to H2-SIM-700084420-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700084420-60 100% 3707KB 741.4KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700084480-60.gwf to H2-SIM-700084480-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700084480-60 100% 3712KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700084540-60.gwf to H2-SIM-700084540-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700084540-60 100% 3716KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700084600-60.gwf to H2-SIM-700084600-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700084600-60 100% 3701KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700084660-60.gwf to H2-SIM-700084660-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700084660-60 100% 3708KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700084720-60.gwf to H2-SIM-700084720-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700084720-60 100% 3715KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700084780-60.gwf to H2-SIM-700084780-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700084780-60 100% 3716KB 371.6KB/s   00:10    
Fetching /home/tania/Popcorn/data2/H2-SIM-700084840-60.gwf to H2-SIM-700084840-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700084840-60 100% 3715KB 128.1KB/s   00:29    
Fetching /home/tania/Popcorn/data2/H2-SIM-700084900-60.gwf to H2-SIM-700084900-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700084900-60 100% 3714KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700084960-60.gwf to H2-SIM-700084960-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700084960-60 100% 3708KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700085020-60.gwf to H2-SIM-700085020-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700085020-60 100% 3712KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700085080-60.gwf to H2-SIM-700085080-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700085080-60 100% 3711KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700085140-60.gwf to H2-SIM-700085140-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700085140-60 100% 3703KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700085200-60.gwf to H2-SIM-700085200-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700085200-60 100% 3700KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700085260-60.gwf to H2-SIM-700085260-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700085260-60 100% 3717KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700085320-60.gwf to H2-SIM-700085320-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700085320-60 100% 3708KB 927.0KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700085380-60.gwf to H2-SIM-700085380-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700085380-60 100% 3716KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700085440-60.gwf to H2-SIM-700085440-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700085440-60 100% 3721KB 930.4KB/s   00:04    
Fetching /home/tania/Popcorn/data2/H2-SIM-700085500-60.gwf to H2-SIM-700085500-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700085500-60 100% 3712KB 742.3KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700085560-60.gwf to H2-SIM-700085560-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700085560-60 100% 3707KB 264.8KB/s   00:14    
Fetching /home/tania/Popcorn/data2/H2-SIM-700085620-60.gwf to H2-SIM-700085620-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700085620-60 100% 3711KB 142.7KB/s   00:26    
Fetching /home/tania/Popcorn/data2/H2-SIM-700085680-60.gwf to H2-SIM-700085680-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700085680-60 100% 3709KB 741.9KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700085740-60.gwf to H2-SIM-700085740-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700085740-60 100% 3718KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700085800-60.gwf to H2-SIM-700085800-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700085800-60 100% 3706KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700085860-60.gwf to H2-SIM-700085860-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700085860-60 100% 3720KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700085920-60.gwf to H2-SIM-700085920-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700085920-60 100% 3711KB   1.8MB/s   00:02    
Fetching /home/tania/Popcorn/data2/H2-SIM-700085980-60.gwf to H2-SIM-700085980-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700085980-60 100% 3712KB   3.6MB/s   00:01    
Fetching /home/tania/Popcorn/data2/H2-SIM-700086040-60.gwf to H2-SIM-700086040-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700086040-60 100% 3704KB   3.6MB/s   00:00    
Fetching /home/tania/Popcorn/data2/H2-SIM-700086100-60.gwf to H2-SIM-700086100-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700086100-60 100% 3713KB 742.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700086160-60.gwf to H2-SIM-700086160-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700086160-60 100% 3706KB   1.2MB/s   00:03    
Fetching /home/tania/Popcorn/data2/H2-SIM-700086220-60.gwf to H2-SIM-700086220-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700086220-60 100% 3713KB 742.5KB/s   00:05    
Fetching /home/tania/Popcorn/data2/H2-SIM-700086280-60.gwf to H2-SIM-700086280-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700086280-60 100% 3713KB 618.8KB/s   00:06    
Fetching /home/tania/Popcorn/data2/H2-SIM-700086340-60.gwf to H2-SIM-700086340-60.gwf
/home/tania/Popcorn/data2/H2-SIM-700086340-60 100% 3719KB 619.9KB/s   00:06    
sftp> sftp> sftp> sftp> sftp> sftp> sftp> sftp> sftp> sftp> sftp> sftp> sftp> sftp> sftp> sftp> sftp> sftp> sftp> sftp> sftp> sftp> exit
[regimbau@projekct ~]$ ls
codes.tar                H2-SIM-700029040-60.gwf  H2-SIM-700058200-60.gwf
data                     H2-SIM-700029100-60.gwf  H2-SIM-700058260-60.gwf
H2-SIM-700000000-60.gwf  H2-SIM-700029160-60.gwf  H2-SIM-700058320-60.gwf
H2-SIM-700000060-60.gwf  H2-SIM-700029220-60.gwf  H2-SIM-700058380-60.gwf
H2-SIM-700000120-60.gwf  H2-SIM-700029280-60.gwf  H2-SIM-700058440-60.gwf
H2-SIM-700000180-60.gwf  H2-SIM-700029340-60.gwf  H2-SIM-700058500-60.gwf
H2-SIM-700000240-60.gwf  H2-SIM-700029400-60.gwf  H2-SIM-700058560-60.gwf
H2-SIM-700000300-60.gwf  H2-SIM-700029460-60.gwf  H2-SIM-700058620-60.gwf
H2-SIM-700000360-60.gwf  H2-SIM-700029520-60.gwf  H2-SIM-700058680-60.gwf
H2-SIM-700000420-60.gwf  H2-SIM-700029580-60.gwf  H2-SIM-700058740-60.gwf
H2-SIM-700000480-60.gwf  H2-SIM-700029640-60.gwf  H2-SIM-700058800-60.gwf
H2-SIM-700000540-60.gwf  H2-SIM-700029700-60.gwf  H2-SIM-700058860-60.gwf
H2-SIM-700000600-60.gwf  H2-SIM-700029760-60.gwf  H2-SIM-700058920-60.gwf
H2-SIM-700000660-60.gwf  H2-SIM-700029820-60.gwf  H2-SIM-700058980-60.gwf
H2-SIM-700000720-60.gwf  H2-SIM-700029880-60.gwf  H2-SIM-700059040-60.gwf
H2-SIM-700000780-60.gwf  H2-SIM-700029940-60.gwf  H2-SIM-700059100-60.gwf
H2-SIM-700000840-60.gwf  H2-SIM-700030000-60.gwf  H2-SIM-700059160-60.gwf
H2-SIM-700000900-60.gwf  H2-SIM-700030060-60.gwf  H2-SIM-700059220-60.gwf
H2-SIM-700000960-60.gwf  H2-SIM-700030120-60.gwf  H2-SIM-700059280-60.gwf
H2-SIM-700001020-60.gwf  H2-SIM-700030180-60.gwf  H2-SIM-700059340-60.gwf
H2-SIM-700001080-60.gwf  H2-SIM-700030240-60.gwf  H2-SIM-700059400-60.gwf
H2-SIM-700001140-60.gwf  H2-SIM-700030300-60.gwf  H2-SIM-700059460-60.gwf
H2-SIM-700001200-60.gwf  H2-SIM-700030360-60.gwf  H2-SIM-700059520-60.gwf
H2-SIM-700001260-60.gwf  H2-SIM-700030420-60.gwf  H2-SIM-700059580-60.gwf
H2-SIM-700001320-60.gwf  H2-SIM-700030480-60.gwf  H2-SIM-700059640-60.gwf
H2-SIM-700001380-60.gwf  H2-SIM-700030540-60.gwf  H2-SIM-700059700-60.gwf
H2-SIM-700001440-60.gwf  H2-SIM-700030600-60.gwf  H2-SIM-700059760-60.gwf
H2-SIM-700001500-60.gwf  H2-SIM-700030660-60.gwf  H2-SIM-700059820-60.gwf
H2-SIM-700001560-60.gwf  H2-SIM-700030720-60.gwf  H2-SIM-700059880-60.gwf
H2-SIM-700001620-60.gwf  H2-SIM-700030780-60.gwf  H2-SIM-700059940-60.gwf
H2-SIM-700001680-60.gwf  H2-SIM-700030840-60.gwf  H2-SIM-700060000-60.gwf
H2-SIM-700001740-60.gwf  H2-SIM-700030900-60.gwf  H2-SIM-700060060-60.gwf
H2-SIM-700001800-60.gwf  H2-SIM-700030960-60.gwf  H2-SIM-700060120-60.gwf
H2-SIM-700001860-60.gwf  H2-SIM-700031020-60.gwf  H2-SIM-700060180-60.gwf
H2-SIM-700001920-60.gwf  H2-SIM-700031080-60.gwf  H2-SIM-700060240-60.gwf
H2-SIM-700001980-60.gwf  H2-SIM-700031140-60.gwf  H2-SIM-700060300-60.gwf
H2-SIM-700002040-60.gwf  H2-SIM-700031200-60.gwf  H2-SIM-700060360-60.gwf
H2-SIM-700002100-60.gwf  H2-SIM-700031260-60.gwf  H2-SIM-700060420-60.gwf
H2-SIM-700002160-60.gwf  H2-SIM-700031320-60.gwf  H2-SIM-700060480-60.gwf
H2-SIM-700002220-60.gwf  H2-SIM-700031380-60.gwf  H2-SIM-700060540-60.gwf
H2-SIM-700002280-60.gwf  H2-SIM-700031440-60.gwf  H2-SIM-700060600-60.gwf
H2-SIM-700002340-60.gwf  H2-SIM-700031500-60.gwf  H2-SIM-700060660-60.gwf
H2-SIM-700002400-60.gwf  H2-SIM-700031560-60.gwf  H2-SIM-700060720-60.gwf
H2-SIM-700002460-60.gwf  H2-SIM-700031620-60.gwf  H2-SIM-700060780-60.gwf
H2-SIM-700002520-60.gwf  H2-SIM-700031680-60.gwf  H2-SIM-700060840-60.gwf
H2-SIM-700002580-60.gwf  H2-SIM-700031740-60.gwf  H2-SIM-700060900-60.gwf
H2-SIM-700002640-60.gwf  H2-SIM-700031800-60.gwf  H2-SIM-700060960-60.gwf
H2-SIM-700002700-60.gwf  H2-SIM-700031860-60.gwf  H2-SIM-700061020-60.gwf
H2-SIM-700002760-60.gwf  H2-SIM-700031920-60.gwf  H2-SIM-700061080-60.gwf
H2-SIM-700002820-60.gwf  H2-SIM-700031980-60.gwf  H2-SIM-700061140-60.gwf
H2-SIM-700002880-60.gwf  H2-SIM-700032040-60.gwf  H2-SIM-700061200-60.gwf
H2-SIM-700002940-60.gwf  H2-SIM-700032100-60.gwf  H2-SIM-700061260-60.gwf
H2-SIM-700003000-60.gwf  H2-SIM-700032160-60.gwf  H2-SIM-700061320-60.gwf
H2-SIM-700003060-60.gwf  H2-SIM-700032220-60.gwf  H2-SIM-700061380-60.gwf
H2-SIM-700003120-60.gwf  H2-SIM-700032280-60.gwf  H2-SIM-700061440-60.gwf
H2-SIM-700003180-60.gwf  H2-SIM-700032340-60.gwf  H2-SIM-700061500-60.gwf
H2-SIM-700003240-60.gwf  H2-SIM-700032400-60.gwf  H2-SIM-700061560-60.gwf
H2-SIM-700003300-60.gwf  H2-SIM-700032460-60.gwf  H2-SIM-700061620-60.gwf
H2-SIM-700003360-60.gwf  H2-SIM-700032520-60.gwf  H2-SIM-700061680-60.gwf
H2-SIM-700003420-60.gwf  H2-SIM-700032580-60.gwf  H2-SIM-700061740-60.gwf
H2-SIM-700003480-60.gwf  H2-SIM-700032640-60.gwf  H2-SIM-700061800-60.gwf
H2-SIM-700003540-60.gwf  H2-SIM-700032700-60.gwf  H2-SIM-700061860-60.gwf
H2-SIM-700003600-60.gwf  H2-SIM-700032760-60.gwf  H2-SIM-700061920-60.gwf
H2-SIM-700003660-60.gwf  H2-SIM-700032820-60.gwf  H2-SIM-700061980-60.gwf
H2-SIM-700003720-60.gwf  H2-SIM-700032880-60.gwf  H2-SIM-700062040-60.gwf
H2-SIM-700003780-60.gwf  H2-SIM-700032940-60.gwf  H2-SIM-700062100-60.gwf
H2-SIM-700003840-60.gwf  H2-SIM-700033000-60.gwf  H2-SIM-700062160-60.gwf
H2-SIM-700003900-60.gwf  H2-SIM-700033060-60.gwf  H2-SIM-700062220-60.gwf
H2-SIM-700003960-60.gwf  H2-SIM-700033120-60.gwf  H2-SIM-700062280-60.gwf
H2-SIM-700004020-60.gwf  H2-SIM-700033180-60.gwf  H2-SIM-700062340-60.gwf
H2-SIM-700004080-60.gwf  H2-SIM-700033240-60.gwf  H2-SIM-700062400-60.gwf
H2-SIM-700004140-60.gwf  H2-SIM-700033300-60.gwf  H2-SIM-700062460-60.gwf
H2-SIM-700004200-60.gwf  H2-SIM-700033360-60.gwf  H2-SIM-700062520-60.gwf
H2-SIM-700004260-60.gwf  H2-SIM-700033420-60.gwf  H2-SIM-700062580-60.gwf
H2-SIM-700004320-60.gwf  H2-SIM-700033480-60.gwf  H2-SIM-700062640-60.gwf
H2-SIM-700004380-60.gwf  H2-SIM-700033540-60.gwf  H2-SIM-700062700-60.gwf
H2-SIM-700004440-60.gwf  H2-SIM-700033600-60.gwf  H2-SIM-700062760-60.gwf
H2-SIM-700004500-60.gwf  H2-SIM-700033660-60.gwf  H2-SIM-700062820-60.gwf
H2-SIM-700004560-60.gwf  H2-SIM-700033720-60.gwf  H2-SIM-700062880-60.gwf
H2-SIM-700004620-60.gwf  H2-SIM-700033780-60.gwf  H2-SIM-700062940-60.gwf
H2-SIM-700004680-60.gwf  H2-SIM-700033840-60.gwf  H2-SIM-700063000-60.gwf
H2-SIM-700004740-60.gwf  H2-SIM-700033900-60.gwf  H2-SIM-700063060-60.gwf
H2-SIM-700004800-60.gwf  H2-SIM-700033960-60.gwf  H2-SIM-700063120-60.gwf
H2-SIM-700004860-60.gwf  H2-SIM-700034020-60.gwf  H2-SIM-700063180-60.gwf
H2-SIM-700004920-60.gwf  H2-SIM-700034080-60.gwf  H2-SIM-700063240-60.gwf
H2-SIM-700004980-60.gwf  H2-SIM-700034140-60.gwf  H2-SIM-700063300-60.gwf
H2-SIM-700005040-60.gwf  H2-SIM-700034200-60.gwf  H2-SIM-700063360-60.gwf
H2-SIM-700005100-60.gwf  H2-SIM-700034260-60.gwf  H2-SIM-700063420-60.gwf
H2-SIM-700005160-60.gwf  H2-SIM-700034320-60.gwf  H2-SIM-700063480-60.gwf
H2-SIM-700005220-60.gwf  H2-SIM-700034380-60.gwf  H2-SIM-700063540-60.gwf
H2-SIM-700005280-60.gwf  H2-SIM-700034440-60.gwf  H2-SIM-700063600-60.gwf
H2-SIM-700005340-60.gwf  H2-SIM-700034500-60.gwf  H2-SIM-700063660-60.gwf
H2-SIM-700005400-60.gwf  H2-SIM-700034560-60.gwf  H2-SIM-700063720-60.gwf
H2-SIM-700005460-60.gwf  H2-SIM-700034620-60.gwf  H2-SIM-700063780-60.gwf
H2-SIM-700005520-60.gwf  H2-SIM-700034680-60.gwf  H2-SIM-700063840-60.gwf
H2-SIM-700005580-60.gwf  H2-SIM-700034740-60.gwf  H2-SIM-700063900-60.gwf
H2-SIM-700005640-60.gwf  H2-SIM-700034800-60.gwf  H2-SIM-700063960-60.gwf
H2-SIM-700005700-60.gwf  H2-SIM-700034860-60.gwf  H2-SIM-700064020-60.gwf
H2-SIM-700005760-60.gwf  H2-SIM-700034920-60.gwf  H2-SIM-700064080-60.gwf
H2-SIM-700005820-60.gwf  H2-SIM-700034980-60.gwf  H2-SIM-700064140-60.gwf
H2-SIM-700005880-60.gwf  H2-SIM-700035040-60.gwf  H2-SIM-700064200-60.gwf
H2-SIM-700005940-60.gwf  H2-SIM-700035100-60.gwf  H2-SIM-700064260-60.gwf
H2-SIM-700006000-60.gwf  H2-SIM-700035160-60.gwf  H2-SIM-700064320-60.gwf
H2-SIM-700006060-60.gwf  H2-SIM-700035220-60.gwf  H2-SIM-700064380-60.gwf
H2-SIM-700006120-60.gwf  H2-SIM-700035280-60.gwf  H2-SIM-700064440-60.gwf
H2-SIM-700006180-60.gwf  H2-SIM-700035340-60.gwf  H2-SIM-700064500-60.gwf
H2-SIM-700006240-60.gwf  H2-SIM-700035400-60.gwf  H2-SIM-700064560-60.gwf
H2-SIM-700006300-60.gwf  H2-SIM-700035460-60.gwf  H2-SIM-700064620-60.gwf
H2-SIM-700006360-60.gwf  H2-SIM-700035520-60.gwf  H2-SIM-700064680-60.gwf
H2-SIM-700006420-60.gwf  H2-SIM-700035580-60.gwf  H2-SIM-700064740-60.gwf
H2-SIM-700006480-60.gwf  H2-SIM-700035640-60.gwf  H2-SIM-700064800-60.gwf
H2-SIM-700006540-60.gwf  H2-SIM-700035700-60.gwf  H2-SIM-700064860-60.gwf
H2-SIM-700006600-60.gwf  H2-SIM-700035760-60.gwf  H2-SIM-700064920-60.gwf
H2-SIM-700006660-60.gwf  H2-SIM-700035820-60.gwf  H2-SIM-700064980-60.gwf
H2-SIM-700006720-60.gwf  H2-SIM-700035880-60.gwf  H2-SIM-700065040-60.gwf
H2-SIM-700006780-60.gwf  H2-SIM-700035940-60.gwf  H2-SIM-700065100-60.gwf
H2-SIM-700006840-60.gwf  H2-SIM-700036000-60.gwf  H2-SIM-700065160-60.gwf
H2-SIM-700006900-60.gwf  H2-SIM-700036060-60.gwf  H2-SIM-700065220-60.gwf
H2-SIM-700006960-60.gwf  H2-SIM-700036120-60.gwf  H2-SIM-700065280-60.gwf
H2-SIM-700007020-60.gwf  H2-SIM-700036180-60.gwf  H2-SIM-700065340-60.gwf
H2-SIM-700007080-60.gwf  H2-SIM-700036240-60.gwf  H2-SIM-700065400-60.gwf
H2-SIM-700007140-60.gwf  H2-SIM-700036300-60.gwf  H2-SIM-700065460-60.gwf
H2-SIM-700007200-60.gwf  H2-SIM-700036360-60.gwf  H2-SIM-700065520-60.gwf
H2-SIM-700007260-60.gwf  H2-SIM-700036420-60.gwf  H2-SIM-700065580-60.gwf
H2-SIM-700007320-60.gwf  H2-SIM-700036480-60.gwf  H2-SIM-700065640-60.gwf
H2-SIM-700007380-60.gwf  H2-SIM-700036540-60.gwf  H2-SIM-700065700-60.gwf
H2-SIM-700007440-60.gwf  H2-SIM-700036600-60.gwf  H2-SIM-700065760-60.gwf
H2-SIM-700007500-60.gwf  H2-SIM-700036660-60.gwf  H2-SIM-700065820-60.gwf
H2-SIM-700007560-60.gwf  H2-SIM-700036720-60.gwf  H2-SIM-700065880-60.gwf
H2-SIM-700007620-60.gwf  H2-SIM-700036780-60.gwf  H2-SIM-700065940-60.gwf
H2-SIM-700007680-60.gwf  H2-SIM-700036840-60.gwf  H2-SIM-700066000-60.gwf
H2-SIM-700007740-60.gwf  H2-SIM-700036900-60.gwf  H2-SIM-700066060-60.gwf
H2-SIM-700007800-60.gwf  H2-SIM-700036960-60.gwf  H2-SIM-700066120-60.gwf
H2-SIM-700007860-60.gwf  H2-SIM-700037020-60.gwf  H2-SIM-700066180-60.gwf
H2-SIM-700007920-60.gwf  H2-SIM-700037080-60.gwf  H2-SIM-700066240-60.gwf
H2-SIM-700007980-60.gwf  H2-SIM-700037140-60.gwf  H2-SIM-700066300-60.gwf
H2-SIM-700008040-60.gwf  H2-SIM-700037200-60.gwf  H2-SIM-700066360-60.gwf
H2-SIM-700008100-60.gwf  H2-SIM-700037260-60.gwf  H2-SIM-700066420-60.gwf
H2-SIM-700008160-60.gwf  H2-SIM-700037320-60.gwf  H2-SIM-700066480-60.gwf
H2-SIM-700008220-60.gwf  H2-SIM-700037380-60.gwf  H2-SIM-700066540-60.gwf
H2-SIM-700008280-60.gwf  H2-SIM-700037440-60.gwf  H2-SIM-700066600-60.gwf
H2-SIM-700008340-60.gwf  H2-SIM-700037500-60.gwf  H2-SIM-700066660-60.gwf
H2-SIM-700008400-60.gwf  H2-SIM-700037560-60.gwf  H2-SIM-700066720-60.gwf
H2-SIM-700008460-60.gwf  H2-SIM-700037620-60.gwf  H2-SIM-700066780-60.gwf
H2-SIM-700008520-60.gwf  H2-SIM-700037680-60.gwf  H2-SIM-700066840-60.gwf
H2-SIM-700008580-60.gwf  H2-SIM-700037740-60.gwf  H2-SIM-700066900-60.gwf
H2-SIM-700008640-60.gwf  H2-SIM-700037800-60.gwf  H2-SIM-700066960-60.gwf
H2-SIM-700008700-60.gwf  H2-SIM-700037860-60.gwf  H2-SIM-700067020-60.gwf
H2-SIM-700008760-60.gwf  H2-SIM-700037920-60.gwf  H2-SIM-700067080-60.gwf
H2-SIM-700008820-60.gwf  H2-SIM-700037980-60.gwf  H2-SIM-700067140-60.gwf
H2-SIM-700008880-60.gwf  H2-SIM-700038040-60.gwf  H2-SIM-700067200-60.gwf
H2-SIM-700008940-60.gwf  H2-SIM-700038100-60.gwf  H2-SIM-700067260-60.gwf
H2-SIM-700009000-60.gwf  H2-SIM-700038160-60.gwf  H2-SIM-700067320-60.gwf
H2-SIM-700009060-60.gwf  H2-SIM-700038220-60.gwf  H2-SIM-700067380-60.gwf
H2-SIM-700009120-60.gwf  H2-SIM-700038280-60.gwf  H2-SIM-700067440-60.gwf
H2-SIM-700009180-60.gwf  H2-SIM-700038340-60.gwf  H2-SIM-700067500-60.gwf
H2-SIM-700009240-60.gwf  H2-SIM-700038400-60.gwf  H2-SIM-700067560-60.gwf
H2-SIM-700009300-60.gwf  H2-SIM-700038460-60.gwf  H2-SIM-700067620-60.gwf
H2-SIM-700009360-60.gwf  H2-SIM-700038520-60.gwf  H2-SIM-700067680-60.gwf
H2-SIM-700009420-60.gwf  H2-SIM-700038580-60.gwf  H2-SIM-700067740-60.gwf
H2-SIM-700009480-60.gwf  H2-SIM-700038640-60.gwf  H2-SIM-700067800-60.gwf
H2-SIM-700009540-60.gwf  H2-SIM-700038700-60.gwf  H2-SIM-700067860-60.gwf
H2-SIM-700009600-60.gwf  H2-SIM-700038760-60.gwf  H2-SIM-700067920-60.gwf
H2-SIM-700009660-60.gwf  H2-SIM-700038820-60.gwf  H2-SIM-700067980-60.gwf
H2-SIM-700009720-60.gwf  H2-SIM-700038880-60.gwf  H2-SIM-700068040-60.gwf
H2-SIM-700009780-60.gwf  H2-SIM-700038940-60.gwf  H2-SIM-700068100-60.gwf
H2-SIM-700009840-60.gwf  H2-SIM-700039000-60.gwf  H2-SIM-700068160-60.gwf
H2-SIM-700009900-60.gwf  H2-SIM-700039060-60.gwf  H2-SIM-700068220-60.gwf
H2-SIM-700009960-60.gwf  H2-SIM-700039120-60.gwf  H2-SIM-700068280-60.gwf
H2-SIM-700010020-60.gwf  H2-SIM-700039180-60.gwf  H2-SIM-700068340-60.gwf
H2-SIM-700010080-60.gwf  H2-SIM-700039240-60.gwf  H2-SIM-700068400-60.gwf
H2-SIM-700010140-60.gwf  H2-SIM-700039300-60.gwf  H2-SIM-700068460-60.gwf
H2-SIM-700010200-60.gwf  H2-SIM-700039360-60.gwf  H2-SIM-700068520-60.gwf
H2-SIM-700010260-60.gwf  H2-SIM-700039420-60.gwf  H2-SIM-700068580-60.gwf
H2-SIM-700010320-60.gwf  H2-SIM-700039480-60.gwf  H2-SIM-700068640-60.gwf
H2-SIM-700010380-60.gwf  H2-SIM-700039540-60.gwf  H2-SIM-700068700-60.gwf
H2-SIM-700010440-60.gwf  H2-SIM-700039600-60.gwf  H2-SIM-700068760-60.gwf
H2-SIM-700010500-60.gwf  H2-SIM-700039660-60.gwf  H2-SIM-700068820-60.gwf
H2-SIM-700010560-60.gwf  H2-SIM-700039720-60.gwf  H2-SIM-700068880-60.gwf
H2-SIM-700010620-60.gwf  H2-SIM-700039780-60.gwf  H2-SIM-700068940-60.gwf
H2-SIM-700010680-60.gwf  H2-SIM-700039840-60.gwf  H2-SIM-700069000-60.gwf
H2-SIM-700010740-60.gwf  H2-SIM-700039900-60.gwf  H2-SIM-700069060-60.gwf
H2-SIM-700010800-60.gwf  H2-SIM-700039960-60.gwf  H2-SIM-700069120-60.gwf
H2-SIM-700010860-60.gwf  H2-SIM-700040020-60.gwf  H2-SIM-700069180-60.gwf
H2-SIM-700010920-60.gwf  H2-SIM-700040080-60.gwf  H2-SIM-700069240-60.gwf
H2-SIM-700010980-60.gwf  H2-SIM-700040140-60.gwf  H2-SIM-700069300-60.gwf
H2-SIM-700011040-60.gwf  H2-SIM-700040200-60.gwf  H2-SIM-700069360-60.gwf
H2-SIM-700011100-60.gwf  H2-SIM-700040260-60.gwf  H2-SIM-700069420-60.gwf
H2-SIM-700011160-60.gwf  H2-SIM-700040320-60.gwf  H2-SIM-700069480-60.gwf
H2-SIM-700011220-60.gwf  H2-SIM-700040380-60.gwf  H2-SIM-700069540-60.gwf
H2-SIM-700011280-60.gwf  H2-SIM-700040440-60.gwf  H2-SIM-700069600-60.gwf
H2-SIM-700011340-60.gwf  H2-SIM-700040500-60.gwf  H2-SIM-700069660-60.gwf
H2-SIM-700011400-60.gwf  H2-SIM-700040560-60.gwf  H2-SIM-700069720-60.gwf
H2-SIM-700011460-60.gwf  H2-SIM-700040620-60.gwf  H2-SIM-700069780-60.gwf
H2-SIM-700011520-60.gwf  H2-SIM-700040680-60.gwf  H2-SIM-700069840-60.gwf
H2-SIM-700011580-60.gwf  H2-SIM-700040740-60.gwf  H2-SIM-700069900-60.gwf
H2-SIM-700011640-60.gwf  H2-SIM-700040800-60.gwf  H2-SIM-700069960-60.gwf
H2-SIM-700011700-60.gwf  H2-SIM-700040860-60.gwf  H2-SIM-700070020-60.gwf
H2-SIM-700011760-60.gwf  H2-SIM-700040920-60.gwf  H2-SIM-700070080-60.gwf
H2-SIM-700011820-60.gwf  H2-SIM-700040980-60.gwf  H2-SIM-700070140-60.gwf
H2-SIM-700011880-60.gwf  H2-SIM-700041040-60.gwf  H2-SIM-700070200-60.gwf
H2-SIM-700011940-60.gwf  H2-SIM-700041100-60.gwf  H2-SIM-700070260-60.gwf
H2-SIM-700012000-60.gwf  H2-SIM-700041160-60.gwf  H2-SIM-700070320-60.gwf
H2-SIM-700012060-60.gwf  H2-SIM-700041220-60.gwf  H2-SIM-700070380-60.gwf
H2-SIM-700012120-60.gwf  H2-SIM-700041280-60.gwf  H2-SIM-700070440-60.gwf
H2-SIM-700012180-60.gwf  H2-SIM-700041340-60.gwf  H2-SIM-700070500-60.gwf
H2-SIM-700012240-60.gwf  H2-SIM-700041400-60.gwf  H2-SIM-700070560-60.gwf
H2-SIM-700012300-60.gwf  H2-SIM-700041460-60.gwf  H2-SIM-700070620-60.gwf
H2-SIM-700012360-60.gwf  H2-SIM-700041520-60.gwf  H2-SIM-700070680-60.gwf
H2-SIM-700012420-60.gwf  H2-SIM-700041580-60.gwf  H2-SIM-700070740-60.gwf
H2-SIM-700012480-60.gwf  H2-SIM-700041640-60.gwf  H2-SIM-700070800-60.gwf
H2-SIM-700012540-60.gwf  H2-SIM-700041700-60.gwf  H2-SIM-700070860-60.gwf
H2-SIM-700012600-60.gwf  H2-SIM-700041760-60.gwf  H2-SIM-700070920-60.gwf
H2-SIM-700012660-60.gwf  H2-SIM-700041820-60.gwf  H2-SIM-700070980-60.gwf
H2-SIM-700012720-60.gwf  H2-SIM-700041880-60.gwf  H2-SIM-700071040-60.gwf
H2-SIM-700012780-60.gwf  H2-SIM-700041940-60.gwf  H2-SIM-700071100-60.gwf
H2-SIM-700012840-60.gwf  H2-SIM-700042000-60.gwf  H2-SIM-700071160-60.gwf
H2-SIM-700012900-60.gwf  H2-SIM-700042060-60.gwf  H2-SIM-700071220-60.gwf
H2-SIM-700012960-60.gwf  H2-SIM-700042120-60.gwf  H2-SIM-700071280-60.gwf
H2-SIM-700013020-60.gwf  H2-SIM-700042180-60.gwf  H2-SIM-700071340-60.gwf
H2-SIM-700013080-60.gwf  H2-SIM-700042240-60.gwf  H2-SIM-700071400-60.gwf
H2-SIM-700013140-60.gwf  H2-SIM-700042300-60.gwf  H2-SIM-700071460-60.gwf
H2-SIM-700013200-60.gwf  H2-SIM-700042360-60.gwf  H2-SIM-700071520-60.gwf
H2-SIM-700013260-60.gwf  H2-SIM-700042420-60.gwf  H2-SIM-700071580-60.gwf
H2-SIM-700013320-60.gwf  H2-SIM-700042480-60.gwf  H2-SIM-700071640-60.gwf
H2-SIM-700013380-60.gwf  H2-SIM-700042540-60.gwf  H2-SIM-700071700-60.gwf
H2-SIM-700013440-60.gwf  H2-SIM-700042600-60.gwf  H2-SIM-700071760-60.gwf
H2-SIM-700013500-60.gwf  H2-SIM-700042660-60.gwf  H2-SIM-700071820-60.gwf
H2-SIM-700013560-60.gwf  H2-SIM-700042720-60.gwf  H2-SIM-700071880-60.gwf
H2-SIM-700013620-60.gwf  H2-SIM-700042780-60.gwf  H2-SIM-700071940-60.gwf
H2-SIM-700013680-60.gwf  H2-SIM-700042840-60.gwf  H2-SIM-700072000-60.gwf
H2-SIM-700013740-60.gwf  H2-SIM-700042900-60.gwf  H2-SIM-700072060-60.gwf
H2-SIM-700013800-60.gwf  H2-SIM-700042960-60.gwf  H2-SIM-700072120-60.gwf
H2-SIM-700013860-60.gwf  H2-SIM-700043020-60.gwf  H2-SIM-700072180-60.gwf
H2-SIM-700013920-60.gwf  H2-SIM-700043080-60.gwf  H2-SIM-700072240-60.gwf
H2-SIM-700013980-60.gwf  H2-SIM-700043140-60.gwf  H2-SIM-700072300-60.gwf
H2-SIM-700014040-60.gwf  H2-SIM-700043200-60.gwf  H2-SIM-700072360-60.gwf
H2-SIM-700014100-60.gwf  H2-SIM-700043260-60.gwf  H2-SIM-700072420-60.gwf
H2-SIM-700014160-60.gwf  H2-SIM-700043320-60.gwf  H2-SIM-700072480-60.gwf
H2-SIM-700014220-60.gwf  H2-SIM-700043380-60.gwf  H2-SIM-700072540-60.gwf
H2-SIM-700014280-60.gwf  H2-SIM-700043440-60.gwf  H2-SIM-700072600-60.gwf
H2-SIM-700014340-60.gwf  H2-SIM-700043500-60.gwf  H2-SIM-700072660-60.gwf
H2-SIM-700014400-60.gwf  H2-SIM-700043560-60.gwf  H2-SIM-700072720-60.gwf
H2-SIM-700014460-60.gwf  H2-SIM-700043620-60.gwf  H2-SIM-700072780-60.gwf
H2-SIM-700014520-60.gwf  H2-SIM-700043680-60.gwf  H2-SIM-700072840-60.gwf
H2-SIM-700014580-60.gwf  H2-SIM-700043740-60.gwf  H2-SIM-700072900-60.gwf
H2-SIM-700014640-60.gwf  H2-SIM-700043800-60.gwf  H2-SIM-700072960-60.gwf
H2-SIM-700014700-60.gwf  H2-SIM-700043860-60.gwf  H2-SIM-700073020-60.gwf
H2-SIM-700014760-60.gwf  H2-SIM-700043920-60.gwf  H2-SIM-700073080-60.gwf
H2-SIM-700014820-60.gwf  H2-SIM-700043980-60.gwf  H2-SIM-700073140-60.gwf
H2-SIM-700014880-60.gwf  H2-SIM-700044040-60.gwf  H2-SIM-700073200-60.gwf
H2-SIM-700014940-60.gwf  H2-SIM-700044100-60.gwf  H2-SIM-700073260-60.gwf
H2-SIM-700015000-60.gwf  H2-SIM-700044160-60.gwf  H2-SIM-700073320-60.gwf
H2-SIM-700015060-60.gwf  H2-SIM-700044220-60.gwf  H2-SIM-700073380-60.gwf
H2-SIM-700015120-60.gwf  H2-SIM-700044280-60.gwf  H2-SIM-700073440-60.gwf
H2-SIM-700015180-60.gwf  H2-SIM-700044340-60.gwf  H2-SIM-700073500-60.gwf
H2-SIM-700015240-60.gwf  H2-SIM-700044400-60.gwf  H2-SIM-700073560-60.gwf
H2-SIM-700015300-60.gwf  H2-SIM-700044460-60.gwf  H2-SIM-700073620-60.gwf
H2-SIM-700015360-60.gwf  H2-SIM-700044520-60.gwf  H2-SIM-700073680-60.gwf
H2-SIM-700015420-60.gwf  H2-SIM-700044580-60.gwf  H2-SIM-700073740-60.gwf
H2-SIM-700015480-60.gwf  H2-SIM-700044640-60.gwf  H2-SIM-700073800-60.gwf
H2-SIM-700015540-60.gwf  H2-SIM-700044700-60.gwf  H2-SIM-700073860-60.gwf
H2-SIM-700015600-60.gwf  H2-SIM-700044760-60.gwf  H2-SIM-700073920-60.gwf
H2-SIM-700015660-60.gwf  H2-SIM-700044820-60.gwf  H2-SIM-700073980-60.gwf
H2-SIM-700015720-60.gwf  H2-SIM-700044880-60.gwf  H2-SIM-700074040-60.gwf
H2-SIM-700015780-60.gwf  H2-SIM-700044940-60.gwf  H2-SIM-700074100-60.gwf
H2-SIM-700015840-60.gwf  H2-SIM-700045000-60.gwf  H2-SIM-700074160-60.gwf
H2-SIM-700015900-60.gwf  H2-SIM-700045060-60.gwf  H2-SIM-700074220-60.gwf
H2-SIM-700015960-60.gwf  H2-SIM-700045120-60.gwf  H2-SIM-700074280-60.gwf
H2-SIM-700016020-60.gwf  H2-SIM-700045180-60.gwf  H2-SIM-700074340-60.gwf
H2-SIM-700016080-60.gwf  H2-SIM-700045240-60.gwf  H2-SIM-700074400-60.gwf
H2-SIM-700016140-60.gwf  H2-SIM-700045300-60.gwf  H2-SIM-700074460-60.gwf
H2-SIM-700016200-60.gwf  H2-SIM-700045360-60.gwf  H2-SIM-700074520-60.gwf
H2-SIM-700016260-60.gwf  H2-SIM-700045420-60.gwf  H2-SIM-700074580-60.gwf
H2-SIM-700016320-60.gwf  H2-SIM-700045480-60.gwf  H2-SIM-700074640-60.gwf
H2-SIM-700016380-60.gwf  H2-SIM-700045540-60.gwf  H2-SIM-700074700-60.gwf
H2-SIM-700016440-60.gwf  H2-SIM-700045600-60.gwf  H2-SIM-700074760-60.gwf
H2-SIM-700016500-60.gwf  H2-SIM-700045660-60.gwf  H2-SIM-700074820-60.gwf
H2-SIM-700016560-60.gwf  H2-SIM-700045720-60.gwf  H2-SIM-700074880-60.gwf
H2-SIM-700016620-60.gwf  H2-SIM-700045780-60.gwf  H2-SIM-700074940-60.gwf
H2-SIM-700016680-60.gwf  H2-SIM-700045840-60.gwf  H2-SIM-700075000-60.gwf
H2-SIM-700016740-60.gwf  H2-SIM-700045900-60.gwf  H2-SIM-700075060-60.gwf
H2-SIM-700016800-60.gwf  H2-SIM-700045960-60.gwf  H2-SIM-700075120-60.gwf
H2-SIM-700016860-60.gwf  H2-SIM-700046020-60.gwf  H2-SIM-700075180-60.gwf
H2-SIM-700016920-60.gwf  H2-SIM-700046080-60.gwf  H2-SIM-700075240-60.gwf
H2-SIM-700016980-60.gwf  H2-SIM-700046140-60.gwf  H2-SIM-700075300-60.gwf
H2-SIM-700017040-60.gwf  H2-SIM-700046200-60.gwf  H2-SIM-700075360-60.gwf
H2-SIM-700017100-60.gwf  H2-SIM-700046260-60.gwf  H2-SIM-700075420-60.gwf
H2-SIM-700017160-60.gwf  H2-SIM-700046320-60.gwf  H2-SIM-700075480-60.gwf
H2-SIM-700017220-60.gwf  H2-SIM-700046380-60.gwf  H2-SIM-700075540-60.gwf
H2-SIM-700017280-60.gwf  H2-SIM-700046440-60.gwf  H2-SIM-700075600-60.gwf
H2-SIM-700017340-60.gwf  H2-SIM-700046500-60.gwf  H2-SIM-700075660-60.gwf
H2-SIM-700017400-60.gwf  H2-SIM-700046560-60.gwf  H2-SIM-700075720-60.gwf
H2-SIM-700017460-60.gwf  H2-SIM-700046620-60.gwf  H2-SIM-700075780-60.gwf
H2-SIM-700017520-60.gwf  H2-SIM-700046680-60.gwf  H2-SIM-700075840-60.gwf
H2-SIM-700017580-60.gwf  H2-SIM-700046740-60.gwf  H2-SIM-700075900-60.gwf
H2-SIM-700017640-60.gwf  H2-SIM-700046800-60.gwf  H2-SIM-700075960-60.gwf
H2-SIM-700017700-60.gwf  H2-SIM-700046860-60.gwf  H2-SIM-700076020-60.gwf
H2-SIM-700017760-60.gwf  H2-SIM-700046920-60.gwf  H2-SIM-700076080-60.gwf
H2-SIM-700017820-60.gwf  H2-SIM-700046980-60.gwf  H2-SIM-700076140-60.gwf
H2-SIM-700017880-60.gwf  H2-SIM-700047040-60.gwf  H2-SIM-700076200-60.gwf
H2-SIM-700017940-60.gwf  H2-SIM-700047100-60.gwf  H2-SIM-700076260-60.gwf
H2-SIM-700018000-60.gwf  H2-SIM-700047160-60.gwf  H2-SIM-700076320-60.gwf
H2-SIM-700018060-60.gwf  H2-SIM-700047220-60.gwf  H2-SIM-700076380-60.gwf
H2-SIM-700018120-60.gwf  H2-SIM-700047280-60.gwf  H2-SIM-700076440-60.gwf
H2-SIM-700018180-60.gwf  H2-SIM-700047340-60.gwf  H2-SIM-700076500-60.gwf
H2-SIM-700018240-60.gwf  H2-SIM-700047400-60.gwf  H2-SIM-700076560-60.gwf
H2-SIM-700018300-60.gwf  H2-SIM-700047460-60.gwf  H2-SIM-700076620-60.gwf
H2-SIM-700018360-60.gwf  H2-SIM-700047520-60.gwf  H2-SIM-700076680-60.gwf
H2-SIM-700018420-60.gwf  H2-SIM-700047580-60.gwf  H2-SIM-700076740-60.gwf
H2-SIM-700018480-60.gwf  H2-SIM-700047640-60.gwf  H2-SIM-700076800-60.gwf
H2-SIM-700018540-60.gwf  H2-SIM-700047700-60.gwf  H2-SIM-700076860-60.gwf
H2-SIM-700018600-60.gwf  H2-SIM-700047760-60.gwf  H2-SIM-700076920-60.gwf
H2-SIM-700018660-60.gwf  H2-SIM-700047820-60.gwf  H2-SIM-700076980-60.gwf
H2-SIM-700018720-60.gwf  H2-SIM-700047880-60.gwf  H2-SIM-700077040-60.gwf
H2-SIM-700018780-60.gwf  H2-SIM-700047940-60.gwf  H2-SIM-700077100-60.gwf
H2-SIM-700018840-60.gwf  H2-SIM-700048000-60.gwf  H2-SIM-700077160-60.gwf
H2-SIM-700018900-60.gwf  H2-SIM-700048060-60.gwf  H2-SIM-700077220-60.gwf
H2-SIM-700018960-60.gwf  H2-SIM-700048120-60.gwf  H2-SIM-700077280-60.gwf
H2-SIM-700019020-60.gwf  H2-SIM-700048180-60.gwf  H2-SIM-700077340-60.gwf
H2-SIM-700019080-60.gwf  H2-SIM-700048240-60.gwf  H2-SIM-700077400-60.gwf
H2-SIM-700019140-60.gwf  H2-SIM-700048300-60.gwf  H2-SIM-700077460-60.gwf
H2-SIM-700019200-60.gwf  H2-SIM-700048360-60.gwf  H2-SIM-700077520-60.gwf
H2-SIM-700019260-60.gwf  H2-SIM-700048420-60.gwf  H2-SIM-700077580-60.gwf
H2-SIM-700019320-60.gwf  H2-SIM-700048480-60.gwf  H2-SIM-700077640-60.gwf
H2-SIM-700019380-60.gwf  H2-SIM-700048540-60.gwf  H2-SIM-700077700-60.gwf
H2-SIM-700019440-60.gwf  H2-SIM-700048600-60.gwf  H2-SIM-700077760-60.gwf
H2-SIM-700019500-60.gwf  H2-SIM-700048660-60.gwf  H2-SIM-700077820-60.gwf
H2-SIM-700019560-60.gwf  H2-SIM-700048720-60.gwf  H2-SIM-700077880-60.gwf
H2-SIM-700019620-60.gwf  H2-SIM-700048780-60.gwf  H2-SIM-700077940-60.gwf
H2-SIM-700019680-60.gwf  H2-SIM-700048840-60.gwf  H2-SIM-700078000-60.gwf
H2-SIM-700019740-60.gwf  H2-SIM-700048900-60.gwf  H2-SIM-700078060-60.gwf
H2-SIM-700019800-60.gwf  H2-SIM-700048960-60.gwf  H2-SIM-700078120-60.gwf
H2-SIM-700019860-60.gwf  H2-SIM-700049020-60.gwf  H2-SIM-700078180-60.gwf
H2-SIM-700019920-60.gwf  H2-SIM-700049080-60.gwf  H2-SIM-700078240-60.gwf
H2-SIM-700019980-60.gwf  H2-SIM-700049140-60.gwf  H2-SIM-700078300-60.gwf
H2-SIM-700020040-60.gwf  H2-SIM-700049200-60.gwf  H2-SIM-700078360-60.gwf
H2-SIM-700020100-60.gwf  H2-SIM-700049260-60.gwf  H2-SIM-700078420-60.gwf
H2-SIM-700020160-60.gwf  H2-SIM-700049320-60.gwf  H2-SIM-700078480-60.gwf
H2-SIM-700020220-60.gwf  H2-SIM-700049380-60.gwf  H2-SIM-700078540-60.gwf
H2-SIM-700020280-60.gwf  H2-SIM-700049440-60.gwf  H2-SIM-700078600-60.gwf
H2-SIM-700020340-60.gwf  H2-SIM-700049500-60.gwf  H2-SIM-700078660-60.gwf
H2-SIM-700020400-60.gwf  H2-SIM-700049560-60.gwf  H2-SIM-700078720-60.gwf
H2-SIM-700020460-60.gwf  H2-SIM-700049620-60.gwf  H2-SIM-700078780-60.gwf
H2-SIM-700020520-60.gwf  H2-SIM-700049680-60.gwf  H2-SIM-700078840-60.gwf
H2-SIM-700020580-60.gwf  H2-SIM-700049740-60.gwf  H2-SIM-700078900-60.gwf
H2-SIM-700020640-60.gwf  H2-SIM-700049800-60.gwf  H2-SIM-700078960-60.gwf
H2-SIM-700020700-60.gwf  H2-SIM-700049860-60.gwf  H2-SIM-700079020-60.gwf
H2-SIM-700020760-60.gwf  H2-SIM-700049920-60.gwf  H2-SIM-700079080-60.gwf
H2-SIM-700020820-60.gwf  H2-SIM-700049980-60.gwf  H2-SIM-700079140-60.gwf
H2-SIM-700020880-60.gwf  H2-SIM-700050040-60.gwf  H2-SIM-700079200-60.gwf
H2-SIM-700020940-60.gwf  H2-SIM-700050100-60.gwf  H2-SIM-700079260-60.gwf
H2-SIM-700021000-60.gwf  H2-SIM-700050160-60.gwf  H2-SIM-700079320-60.gwf
H2-SIM-700021060-60.gwf  H2-SIM-700050220-60.gwf  H2-SIM-700079380-60.gwf
H2-SIM-700021120-60.gwf  H2-SIM-700050280-60.gwf  H2-SIM-700079440-60.gwf
H2-SIM-700021180-60.gwf  H2-SIM-700050340-60.gwf  H2-SIM-700079500-60.gwf
H2-SIM-700021240-60.gwf  H2-SIM-700050400-60.gwf  H2-SIM-700079560-60.gwf
H2-SIM-700021300-60.gwf  H2-SIM-700050460-60.gwf  H2-SIM-700079620-60.gwf
H2-SIM-700021360-60.gwf  H2-SIM-700050520-60.gwf  H2-SIM-700079680-60.gwf
H2-SIM-700021420-60.gwf  H2-SIM-700050580-60.gwf  H2-SIM-700079740-60.gwf
H2-SIM-700021480-60.gwf  H2-SIM-700050640-60.gwf  H2-SIM-700079800-60.gwf
H2-SIM-700021540-60.gwf  H2-SIM-700050700-60.gwf  H2-SIM-700079860-60.gwf
H2-SIM-700021600-60.gwf  H2-SIM-700050760-60.gwf  H2-SIM-700079920-60.gwf
H2-SIM-700021660-60.gwf  H2-SIM-700050820-60.gwf  H2-SIM-700079980-60.gwf
H2-SIM-700021720-60.gwf  H2-SIM-700050880-60.gwf  H2-SIM-700080040-60.gwf
H2-SIM-700021780-60.gwf  H2-SIM-700050940-60.gwf  H2-SIM-700080100-60.gwf
H2-SIM-700021840-60.gwf  H2-SIM-700051000-60.gwf  H2-SIM-700080160-60.gwf
H2-SIM-700021900-60.gwf  H2-SIM-700051060-60.gwf  H2-SIM-700080220-60.gwf
H2-SIM-700021960-60.gwf  H2-SIM-700051120-60.gwf  H2-SIM-700080280-60.gwf
H2-SIM-700022020-60.gwf  H2-SIM-700051180-60.gwf  H2-SIM-700080340-60.gwf
H2-SIM-700022080-60.gwf  H2-SIM-700051240-60.gwf  H2-SIM-700080400-60.gwf
H2-SIM-700022140-60.gwf  H2-SIM-700051300-60.gwf  H2-SIM-700080460-60.gwf
H2-SIM-700022200-60.gwf  H2-SIM-700051360-60.gwf  H2-SIM-700080520-60.gwf
H2-SIM-700022260-60.gwf  H2-SIM-700051420-60.gwf  H2-SIM-700080580-60.gwf
H2-SIM-700022320-60.gwf  H2-SIM-700051480-60.gwf  H2-SIM-700080640-60.gwf
H2-SIM-700022380-60.gwf  H2-SIM-700051540-60.gwf  H2-SIM-700080700-60.gwf
H2-SIM-700022440-60.gwf  H2-SIM-700051600-60.gwf  H2-SIM-700080760-60.gwf
H2-SIM-700022500-60.gwf  H2-SIM-700051660-60.gwf  H2-SIM-700080820-60.gwf
H2-SIM-700022560-60.gwf  H2-SIM-700051720-60.gwf  H2-SIM-700080880-60.gwf
H2-SIM-700022620-60.gwf  H2-SIM-700051780-60.gwf  H2-SIM-700080940-60.gwf
H2-SIM-700022680-60.gwf  H2-SIM-700051840-60.gwf  H2-SIM-700081000-60.gwf
H2-SIM-700022740-60.gwf  H2-SIM-700051900-60.gwf  H2-SIM-700081060-60.gwf
H2-SIM-700022800-60.gwf  H2-SIM-700051960-60.gwf  H2-SIM-700081120-60.gwf
H2-SIM-700022860-60.gwf  H2-SIM-700052020-60.gwf  H2-SIM-700081180-60.gwf
H2-SIM-700022920-60.gwf  H2-SIM-700052080-60.gwf  H2-SIM-700081240-60.gwf
H2-SIM-700022980-60.gwf  H2-SIM-700052140-60.gwf  H2-SIM-700081300-60.gwf
H2-SIM-700023040-60.gwf  H2-SIM-700052200-60.gwf  H2-SIM-700081360-60.gwf
H2-SIM-700023100-60.gwf  H2-SIM-700052260-60.gwf  H2-SIM-700081420-60.gwf
H2-SIM-700023160-60.gwf  H2-SIM-700052320-60.gwf  H2-SIM-700081480-60.gwf
H2-SIM-700023220-60.gwf  H2-SIM-700052380-60.gwf  H2-SIM-700081540-60.gwf
H2-SIM-700023280-60.gwf  H2-SIM-700052440-60.gwf  H2-SIM-700081600-60.gwf
H2-SIM-700023340-60.gwf  H2-SIM-700052500-60.gwf  H2-SIM-700081660-60.gwf
H2-SIM-700023400-60.gwf  H2-SIM-700052560-60.gwf  H2-SIM-700081720-60.gwf
H2-SIM-700023460-60.gwf  H2-SIM-700052620-60.gwf  H2-SIM-700081780-60.gwf
H2-SIM-700023520-60.gwf  H2-SIM-700052680-60.gwf  H2-SIM-700081840-60.gwf
H2-SIM-700023580-60.gwf  H2-SIM-700052740-60.gwf  H2-SIM-700081900-60.gwf
H2-SIM-700023640-60.gwf  H2-SIM-700052800-60.gwf  H2-SIM-700081960-60.gwf
H2-SIM-700023700-60.gwf  H2-SIM-700052860-60.gwf  H2-SIM-700082020-60.gwf
H2-SIM-700023760-60.gwf  H2-SIM-700052920-60.gwf  H2-SIM-700082080-60.gwf
H2-SIM-700023820-60.gwf  H2-SIM-700052980-60.gwf  H2-SIM-700082140-60.gwf
H2-SIM-700023880-60.gwf  H2-SIM-700053040-60.gwf  H2-SIM-700082200-60.gwf
H2-SIM-700023940-60.gwf  H2-SIM-700053100-60.gwf  H2-SIM-700082260-60.gwf
H2-SIM-700024000-60.gwf  H2-SIM-700053160-60.gwf  H2-SIM-700082320-60.gwf
H2-SIM-700024060-60.gwf  H2-SIM-700053220-60.gwf  H2-SIM-700082380-60.gwf
H2-SIM-700024120-60.gwf  H2-SIM-700053280-60.gwf  H2-SIM-700082440-60.gwf
H2-SIM-700024180-60.gwf  H2-SIM-700053340-60.gwf  H2-SIM-700082500-60.gwf
H2-SIM-700024240-60.gwf  H2-SIM-700053400-60.gwf  H2-SIM-700082560-60.gwf
H2-SIM-700024300-60.gwf  H2-SIM-700053460-60.gwf  H2-SIM-700082620-60.gwf
H2-SIM-700024360-60.gwf  H2-SIM-700053520-60.gwf  H2-SIM-700082680-60.gwf
H2-SIM-700024420-60.gwf  H2-SIM-700053580-60.gwf  H2-SIM-700082740-60.gwf
H2-SIM-700024480-60.gwf  H2-SIM-700053640-60.gwf  H2-SIM-700082800-60.gwf
H2-SIM-700024540-60.gwf  H2-SIM-700053700-60.gwf  H2-SIM-700082860-60.gwf
H2-SIM-700024600-60.gwf  H2-SIM-700053760-60.gwf  H2-SIM-700082920-60.gwf
H2-SIM-700024660-60.gwf  H2-SIM-700053820-60.gwf  H2-SIM-700082980-60.gwf
H2-SIM-700024720-60.gwf  H2-SIM-700053880-60.gwf  H2-SIM-700083040-60.gwf
H2-SIM-700024780-60.gwf  H2-SIM-700053940-60.gwf  H2-SIM-700083100-60.gwf
H2-SIM-700024840-60.gwf  H2-SIM-700054000-60.gwf  H2-SIM-700083160-60.gwf
H2-SIM-700024900-60.gwf  H2-SIM-700054060-60.gwf  H2-SIM-700083220-60.gwf
H2-SIM-700024960-60.gwf  H2-SIM-700054120-60.gwf  H2-SIM-700083280-60.gwf
H2-SIM-700025020-60.gwf  H2-SIM-700054180-60.gwf  H2-SIM-700083340-60.gwf
H2-SIM-700025080-60.gwf  H2-SIM-700054240-60.gwf  H2-SIM-700083400-60.gwf
H2-SIM-700025140-60.gwf  H2-SIM-700054300-60.gwf  H2-SIM-700083460-60.gwf
H2-SIM-700025200-60.gwf  H2-SIM-700054360-60.gwf  H2-SIM-700083520-60.gwf
H2-SIM-700025260-60.gwf  H2-SIM-700054420-60.gwf  H2-SIM-700083580-60.gwf
H2-SIM-700025320-60.gwf  H2-SIM-700054480-60.gwf  H2-SIM-700083640-60.gwf
H2-SIM-700025380-60.gwf  H2-SIM-700054540-60.gwf  H2-SIM-700083700-60.gwf
H2-SIM-700025440-60.gwf  H2-SIM-700054600-60.gwf  H2-SIM-700083760-60.gwf
H2-SIM-700025500-60.gwf  H2-SIM-700054660-60.gwf  H2-SIM-700083820-60.gwf
H2-SIM-700025560-60.gwf  H2-SIM-700054720-60.gwf  H2-SIM-700083880-60.gwf
H2-SIM-700025620-60.gwf  H2-SIM-700054780-60.gwf  H2-SIM-700083940-60.gwf
H2-SIM-700025680-60.gwf  H2-SIM-700054840-60.gwf  H2-SIM-700084000-60.gwf
H2-SIM-700025740-60.gwf  H2-SIM-700054900-60.gwf  H2-SIM-700084060-60.gwf
H2-SIM-700025800-60.gwf  H2-SIM-700054960-60.gwf  H2-SIM-700084120-60.gwf
H2-SIM-700025860-60.gwf  H2-SIM-700055020-60.gwf  H2-SIM-700084180-60.gwf
H2-SIM-700025920-60.gwf  H2-SIM-700055080-60.gwf  H2-SIM-700084240-60.gwf
H2-SIM-700025980-60.gwf  H2-SIM-700055140-60.gwf  H2-SIM-700084300-60.gwf
H2-SIM-700026040-60.gwf  H2-SIM-700055200-60.gwf  H2-SIM-700084360-60.gwf
H2-SIM-700026100-60.gwf  H2-SIM-700055260-60.gwf  H2-SIM-700084420-60.gwf
H2-SIM-700026160-60.gwf  H2-SIM-700055320-60.gwf  H2-SIM-700084480-60.gwf
H2-SIM-700026220-60.gwf  H2-SIM-700055380-60.gwf  H2-SIM-700084540-60.gwf
H2-SIM-700026280-60.gwf  H2-SIM-700055440-60.gwf  H2-SIM-700084600-60.gwf
H2-SIM-700026340-60.gwf  H2-SIM-700055500-60.gwf  H2-SIM-700084660-60.gwf
H2-SIM-700026400-60.gwf  H2-SIM-700055560-60.gwf  H2-SIM-700084720-60.gwf
H2-SIM-700026460-60.gwf  H2-SIM-700055620-60.gwf  H2-SIM-700084780-60.gwf
H2-SIM-700026520-60.gwf  H2-SIM-700055680-60.gwf  H2-SIM-700084840-60.gwf
H2-SIM-700026580-60.gwf  H2-SIM-700055740-60.gwf  H2-SIM-700084900-60.gwf
H2-SIM-700026640-60.gwf  H2-SIM-700055800-60.gwf  H2-SIM-700084960-60.gwf
H2-SIM-700026700-60.gwf  H2-SIM-700055860-60.gwf  H2-SIM-700085020-60.gwf
H2-SIM-700026760-60.gwf  H2-SIM-700055920-60.gwf  H2-SIM-700085080-60.gwf
H2-SIM-700026820-60.gwf  H2-SIM-700055980-60.gwf  H2-SIM-700085140-60.gwf
H2-SIM-700026880-60.gwf  H2-SIM-700056040-60.gwf  H2-SIM-700085200-60.gwf
H2-SIM-700026940-60.gwf  H2-SIM-700056100-60.gwf  H2-SIM-700085260-60.gwf
H2-SIM-700027000-60.gwf  H2-SIM-700056160-60.gwf  H2-SIM-700085320-60.gwf
H2-SIM-700027060-60.gwf  H2-SIM-700056220-60.gwf  H2-SIM-700085380-60.gwf
H2-SIM-700027120-60.gwf  H2-SIM-700056280-60.gwf  H2-SIM-700085440-60.gwf
H2-SIM-700027180-60.gwf  H2-SIM-700056340-60.gwf  H2-SIM-700085500-60.gwf
H2-SIM-700027240-60.gwf  H2-SIM-700056400-60.gwf  H2-SIM-700085560-60.gwf
H2-SIM-700027300-60.gwf  H2-SIM-700056460-60.gwf  H2-SIM-700085620-60.gwf
H2-SIM-700027360-60.gwf  H2-SIM-700056520-60.gwf  H2-SIM-700085680-60.gwf
H2-SIM-700027420-60.gwf  H2-SIM-700056580-60.gwf  H2-SIM-700085740-60.gwf
H2-SIM-700027480-60.gwf  H2-SIM-700056640-60.gwf  H2-SIM-700085800-60.gwf
H2-SIM-700027540-60.gwf  H2-SIM-700056700-60.gwf  H2-SIM-700085860-60.gwf
H2-SIM-700027600-60.gwf  H2-SIM-700056760-60.gwf  H2-SIM-700085920-60.gwf
H2-SIM-700027660-60.gwf  H2-SIM-700056820-60.gwf  H2-SIM-700085980-60.gwf
H2-SIM-700027720-60.gwf  H2-SIM-700056880-60.gwf  H2-SIM-700086040-60.gwf
H2-SIM-700027780-60.gwf  H2-SIM-700056940-60.gwf  H2-SIM-700086100-60.gwf
H2-SIM-700027840-60.gwf  H2-SIM-700057000-60.gwf  H2-SIM-700086160-60.gwf
H2-SIM-700027900-60.gwf  H2-SIM-700057060-60.gwf  H2-SIM-700086220-60.gwf
H2-SIM-700027960-60.gwf  H2-SIM-700057120-60.gwf  H2-SIM-700086280-60.gwf
H2-SIM-700028020-60.gwf  H2-SIM-700057180-60.gwf  H2-SIM-700086340-60.gwf
H2-SIM-700028080-60.gwf  H2-SIM-700057240-60.gwf  ldg-4.4
H2-SIM-700028140-60.gwf  H2-SIM-700057300-60.gwf  nova.lsf~
H2-SIM-700028200-60.gwf  H2-SIM-700057360-60.gwf  opt
H2-SIM-700028260-60.gwf  H2-SIM-700057420-60.gwf  opt.tar
H2-SIM-700028320-60.gwf  H2-SIM-700057480-60.gwf  pacman-3.19
H2-SIM-700028380-60.gwf  H2-SIM-700057540-60.gwf  Popcorn
H2-SIM-700028440-60.gwf  H2-SIM-700057600-60.gwf  popcorn.lsf~
H2-SIM-700028500-60.gwf  H2-SIM-700057660-60.gwf  Popcorn_restore
H2-SIM-700028560-60.gwf  H2-SIM-700057720-60.gwf  pulsar
H2-SIM-700028620-60.gwf  H2-SIM-700057780-60.gwf  regimbaucer.pem
H2-SIM-700028680-60.gwf  H2-SIM-700057840-60.gwf  regimbaukey.pem
H2-SIM-700028740-60.gwf  H2-SIM-700057900-60.gwf  src
H2-SIM-700028800-60.gwf  H2-SIM-700057960-60.gwf  src.tar
H2-SIM-700028860-60.gwf  H2-SIM-700058020-60.gwf  Strings.tar
H2-SIM-700028920-60.gwf  H2-SIM-700058080-60.gwf  supernova
H2-SIM-700028980-60.gwf  H2-SIM-700058140-60.gwf
[regimbau@projekct ~]$ mkdir data2
[regimbau@projekct ~]$ mv H2* data2
[regimbau@projekct ~]$ mv data2 Popcorn
[regimbau@projekct ~]$ cd Popcorn
[regimbau@projekct Popcorn]$ ls
cache.c  gaussian.lsf   popcorn          popcorn.37576.37576  popcorn.c_save
compil   gaussian.lsf~  popcorn.1.37645  popcorn.37611.37611  popcorn.lsf
compil~  H1.cache       popcorn.1.37754  popcorn.37638.1      popcorn.lsf~
compil2  H2.cache       popcorn.2.37645  popcorn.37638.2      popcorn.lsf_save
data1    matlab         popcorn.2.37754  popcorn.c
data2    output         popcorn.3.37754  popcorn.c~
[regimbau@projekct Popcorn]$ bsub -J "popcorn[1-3]" < popcorn.lsfJob <37769> is submitted to default queue <normal>.[regimbau@projekct Popcorn]$ bjobs
JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
37769   regimba RUN   normal     lsfhost.loc lsfhost.loc popcorn[1] Jul 11 23:34
37769   regimba RUN   normal     lsfhost.loc lsfhost.loc popcorn[2] Jul 11 23:34
37769   regimba RUN   normal     lsfhost.loc lsfhost.loc popcorn[3] Jul 11 23:34
[regimbau@projekct Popcorn]$ bsub -J "popcorn[1-3]" < popcorn.lsf
^[[BJob <37770> is submitted to default queue <normal>.
[regimbau@projekct Popcorn]$ bjobs
JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
37769   regimba RUN   normal     lsfhost.loc lsfhost.loc popcorn[1] Jul 11 23:34
37769   regimba RUN   normal     lsfhost.loc lsfhost.loc popcorn[2] Jul 11 23:34
37769   regimba RUN   normal     lsfhost.loc lsfhost.loc popcorn[3] Jul 11 23:34
37770   regimba RUN   normal     lsfhost.loc lsfhost.loc popcorn[1] Jul 11 23:39
37770   regimba RUN   normal     lsfhost.loc lsfhost.loc popcorn[2] Jul 11 23:39
37770   regimba RUN   normal     lsfhost.loc lsfhost.loc popcorn[3] Jul 11 23:39
[regimbau@projekct Popcorn]$ bkill 37770
Job <37770>: Operation is in progress
[regimbau@projekct Popcorn]$ ls output/
stat_1.dat  stat_2.dat  stat_3.dat
[regimbau@projekct Popcorn]$ bjobs
JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
37769   regimba RUN   normal     lsfhost.loc lsfhost.loc popcorn[1] Jul 11 23:34
[regimbau@projekct Popcorn]$ cd output
[regimbau@projekct output]$ more stat_2.dat 
700000300 0.718052 1.392657 0.522510 0.488199 0.124534 0.718128 1.392751
700000361 0.718710 1.391382 1.035171 0.376438 0.146689 0.703291 1.383161
700000422 0.708189 1.412051 0.433681 0.558140 0.135100 0.702075 1.408414
700000483 0.727823 1.373960 0.633047 0.469795 0.139718 0.718299 1.367944
File Edit Options Buffers Tools C Help                                          
/*                                                                              
*  Copyright (C) 2007 Tania Regimbau                                            
*                                                                               
*  This program is free software; you can redistribute it and/or modify         
*  it under the terms of the GNU General Public License as published by         
*  the Free Software Foundation; either version 2 of the License, or            
*  (at your option) any later version.                                          
*                                                                               
*  This program is distributed in the hope that it will be useful,              
*  but WITHOUT ANY WARRANTY; without even the implied warranty of               
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
*  GNU General Public License for more details.                                 
*                                                                               
*  You should have received a copy of the GNU General Public License            
*  along with with program; see the file COPYING. If not, write to the          
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,          
*  MA  02111-1307  USA                                                          
*/

/*                                                                              
 * popcorn.c                                                                    
-uu-:---F1  popcorn.c         (C Abbrev)--L1--Top-------------------------------
Loading cc-mode...done
