"""
followup Web page classes

This
"""

__author__ = 'Chad Hanna <channa@phys.lsu.edu>'

import ConfigParser
import os,re # added for cacheStructure

from pylal import git_version

##############################################################################
# Cachefying class
##############################################################################

class cacheStructure(object):

  def __init__(self):
    self.mainCache = {}

  def appendCache(self,cacheType,outputPath=None):
    try:
      self.mainCache.keys().index(cacheType)
    except:
      self.mainCache[cacheType] = [[],outputPath]
      self.mainCache[cacheType][0].append('### ' + cacheType + ' CACHE FILE ###\n')

  def appendSubCache(self,cacheType,cacheEntry):
    self.mainCache[cacheType][0].append(cacheEntry)

  def writeCache(self):
    if len(self.mainCache) > 0:
      mainFile = open('main.cache','w')
      for subcache in self.mainCache.items():
        fileName = self.fillMainCache(subcache,mainFile)
        if len(subcache[1][0]) > 1:
          cache_file = open(fileName,'w')
          for cacheEntry in subcache[1][0]:
            cache_file.write(cacheEntry)
          cache_file.close()
      mainFile.close()
    else: pass

  def fillMainCache(self,subcache,fileToWrite):
    fileNames = []
    if subcache[1][1]:
      fileName = subcache[1][1] + subcache[0] + '.cache'
    else:
      fileName = subcache[0] + '.cache'
    if subcache[1][1]:
      self.getCacheList(fileNames,subcache[1][1])
    try:
      fileNames.index(fileName)
    except:
      if len(subcache[1][0]) > 1:
        fileNames.append(fileName)
    for file in fileNames:
      fileToWrite.write(subcache[0] + ' ' + file + '\n') 
    return fileName

  def getCacheList(self,list,inputdir):
    list_files = os.listdir(inputdir)
    for element in list_files:
      if re.search('.cache',element):
        list.append(inputdir + element)
      else: pass

##############################################################################
# The webifying stuff starts here
##############################################################################

class Content:

  def __init__(self):
    self.contentList = []
    self.table = []
    self.lastTable = None
    #self.root = ''

  def link(self,link,text):
    thisLink = Link(link,text)
    self.contentList.append(thisLink)

  def verbatim(self,text):
    thisVerbatim = Verbatim(text)
    self.contentList.append(thisVerbatim)

  def text(self,text,type='',color=''):
    thisText = Text(text,type,color)
    self.contentList.append(thisText)

  def list(self,list):
    thisList = List(list)
    self.contentList.append(thisList)

  def image(self, image, link=None, width = 400):
    thisImage = Image(image, self.root, link, width)
    self.contentList.append(thisImage)

  def linebreak(self, times = 1):
    thisBreak = Break(times)
    self.contentList.append(thisBreak)

  def appendTable(self,rows,columns,border=0,width=800):
    number = len(self.table)
    thisTable = Table(rows,columns,border,width,number,self.root)
    self.contentList.append(thisTable)
    self.table.append(thisTable)
    self.lastTable  = self.table[number]


# Currently the write/read methods for talkBack don't work for lists or tables.
class talkBack(Content):

  def __init__(self,outputWebFile):
    self.fileName = outputWebFile.replace('html','ini')
    self.summaryPlot = ''
    self.summaryPlotCaption = ''
    self.summaryText = ''
    self.contentList = []
  def addSummaryPlot(self,plot,caption):
    self.summaryPlot = plot
    self.summaryPlotCaption = caption
  def addSummaryText(self,text):
    self.summaryText = text

  def write(self):
    file = open(self.fileName,'w')
    file.write('[talkBack]\n')
    file.write('\nsummaryPlot='+self.summaryPlot)
    file.write('\nsummaryPlotCaption='+self.summaryPlotCaption)
    file.write('\nsummaryText='+self.summaryText+'\n')
    for content in range(len(self.contentList)):
      self.contentList[content].write(file, 'talkBack',content)
    file.close()
 
  def read(self):
    cp = ConfigParser.ConfigParser()
    cp.read(self.fileName)
    try: self.summaryPlot = cp.get('talkBack','summaryPlot')
    except: pass
    try: self.summaryPlotCaption = cp.get('talkBack','summaryPlotCaption')
    except: pass
    try: self.summaryText = cp.get('talkBack','summaryText')
    except: pass

  #append talkback stuff to a valid webpage class (or section/subsection/table)
  def readAppend(self,webpage):
    cp = ConfigParser.ConfigParser()
    cp.read(self.fileName)
    for section in cp.sections():
      if section.endswith('Link'):
        try: link = cp.get(section,'link')
        except: pass
        try: text = cp.get(section,'text')
        except: pass
        webpage.link(link,text)
      if section.endswith('Text'):
        type = None
        color = None
        try: text = cp.get(section,'text')
        except: pass
        try: type = cp.get(section,'type')
        except: pass
        try: color = cp.get(section,'color')
        except: pass
        webpage.text(text,type,color)
      if section.endswith('Verbatim'):
        try: text = cp.get(section,'text')
        except: pass
        webpage.verbatim(text)
      if section.endswith('Image'):
        try: link = cp.get(section,'link')
        except: link=None
        try: image = cp.get(section,'image')
        except: pass
        try: width = cp.get(section,'width')
        except: width = 400
        webpage.image(image,link,width)
      if section.endswith('Linebreak'):
        try: times = cp.get(section,'linebreak')
        except: pass
        webpage.linebreak(times)


#Add sub pages and provide slick mechanism for writing the whole tree!?!
#this will help with the followup dag web generation...
# ADD CONCEPT OF RELATIVE PATH TO GIVE MORE FLEXIBILITY TO DIRECTORY STRUCTURE
class WebPage(Content):
  """ 
  Class to store the web page output of a followup program
  """
  def __init__(self,title,filename,root=''):
    #Content will be written before sections, which themselves may have content
    self.contentList = []
    self.table = []
    self.section = []
    self.lastSection = None
    self.title = title
    self.root = root
    self.subPage = []
    self.lastPage = None
    self.filename = filename

  def appendSection(self,heading):
    number = len(self.section)
    self.section.append( Section( heading,number,self.filename,self.root ) ) 
    self.lastSection = self.section[number] 

  
  def appendSubPage(self,title, file, root=''):
    self.subPage.append(WebPage(title, file, root))
    self.lastPage = self.subPage[len(self.subPage)]

  def linkNewSubPage(self,title,file, text='', root=''):
    if text:
      self.link(root+file,text)
    else:
      self.link(root+file, title)
    self.appendSubPage(title,file, root)


  def write(self,type):
    self.writeHeader(self.file,type)
    self.writeTitle(self.file,type)
    self.writeTableOfContents(self.file,type)
    # first do the content in this page
    for content in self.contentList:
      content.write(self.file,type)
    for section in self.section:
      section.write(self.file,type)
    # now do the sub pages recursively
    for page in self.subPage:
      page.cleanWrite(type)
    

  def writeTitle(self, file, type):
    if type == 'IUL':
      pass # cause this is done in the header for IUL type

  def writeHeader(self, file, type):
    if type == 'IUL':
      file.write('<%method title>' + self.title + '</%method>\n')
      file.write('<%method headline>' + self.title + '</%method>\n')
      file.write('<%method cvsid>' + git_version.id + '</%method>\n')
      #file.write('<h1>'+self.title+'</h1>\n')
 
  def writeTableOfContents(self,file,type):
    if type == 'IUL':
      file.write('<h2 id="fuwebtoc">Table of contents</h2>\n') 
      sectionTOC  = [] 

      for section in self.section:
        link = section.toclink
        subSectionTOC = []
        for subsection in section.subSection:
          subSectionTOC.append( [Link(subsection.toclink, subsection.heading)])
        if subSectionTOC:
          sectionTOC.append( [Link(section.toclink, section.heading), List(subSectionTOC)] )
        else: 
          sectionTOC.append( [Link(section.toclink, section.heading)] )
      TOC = List(sectionTOC)
      TOC.write(file,type)
          
#  def cleanWrite(self, filename, type):
  def cleanWrite(self,type):
    self.file = open(self.filename,'w')
    self.write(type)          
    self.file.close()

# This class shouldn't really be used without a webpage as it's parent
class Section(Content):
  """
  Class to store a section of a webpage
  """
  def __init__(self,heading,secNumber,filename,root=''):
    self.contentList = []
    self.table = []
    self.subSection = []
    self.lastSub = None
    self.heading = heading
    self.secNumber = secNumber 
    #self.toclink = root + '/'+filename.split('/')[-1]+'#section' + str(self.secNumber)
    self.toclink = root + filename.split('/')[-1]+'#section' + str(self.secNumber)
    self.root = root
    self.filename  = filename

  def appendSubSection(self,heading):
    number = len(self.subSection)
    self.subSection.append( SubSection( heading,self.secNumber,number, self.filename,self.root ) )
    self.lastSub = self.subSection[number]

  def write(self,file,type):
    self.writeSectionHeader(file,type)
    for content in self.contentList:
      content.write(file,type)
    for subSection in self.subSection:
      subSection.write(file,type)
  
  def writeSectionHeader(self,file,type):
    if type == 'IUL':
      file.write('<h2 id="section'+str(self.secNumber)+'">'+str(self.secNumber+1)+'.  ' + self.heading+'\n')
      #file.write('<a href="'+self.root+'/'+file.name.split('/')[-1]+'#fuwebtoc">[Back to TOC]</a></h2>\n')
      file.write('<a href="'+self.root+file.name.split('/')[-1]+'#fuwebtoc">[Back to TOC]</a></h2>\n')
      
# This class shouldn't really be used without a section as its parent, which
# itself has a webpage as its parent
class SubSection(Content):
  """
  Class to store a subsection of a webpage
  """
  def __init__(self,heading,secNumber,subNumber, filename,root=''):
    self.contentList = []
    self.table = []
    self.heading = heading
    self.secNumber = secNumber
    self.subNumber = subNumber
    self.root = root
    self.filename = filename
    #self.toclink = root + '/'+filename.split('/')[-1]+'#subsection' + str(self.secNumber) + '.' + str(self.subNumber)
    self.toclink = root + filename.split('/')[-1]+'#subsection' + str(self.secNumber) + '.' + str(self.subNumber)

  def write(self,file,type):
    self.writeSubSectionHeader(file,type)
    for content in self.contentList:
      content.write(file,type)

  def writeSubSectionHeader(self,file,type):
    if type == 'IUL':
      file.write('<h3 id="subsection'+str(self.secNumber)+'.'+str(self.subNumber)+'">'+str(self.secNumber+1)+'.'+str(self.subNumber+1)+'.  '+self.heading+'\n')
      #file.write('<a href="'+self.root+'/'+file.name.split('/')[-1]+'#fuwebtoc">[Back to TOC]</a></h3>\n')
      file.write('<a href="'+self.root+file.name.split('/')[-1]+'#fuwebtoc">[Back to TOC]</a></h3>\n')



# here is where the real 'content' is.  Remember that all content is
# contained inside of a table. Table allows you to create rows, columns
# and cells, there shouldn't be any normal reason to use those classes
# by themselves  - currently this just does tables of form (m x n) 
class Table:
  """
  Class to store the web page output of a followup program
  """
  def __init__(self,rows,columns,border,width,number,root):
    self.row = []
    self.number = number
    self.border = border
    self.width = width
    self.root = root
    for i in range(rows):
      self.row.append(Row(columns,root))

  def write(self,file,type):
    if type == 'IUL':
      file.write('<table border=' + str(self.border)+' width='+str(self.width)+'>\n')
      for row in self.row:
        row.write(file,type)
      file.write('</table>\n')

# You shouldn't need to use this class - best to use Table
class Row:
  """
  Class to store a table row
  """
  def __init__(self,columns,root):
    self.cell = []
    self.root = root
    for j in range(columns):
      self.cell.append(Cell(root))

  def write(self,file,type):
    if type == 'IUL':
      file.write('<tr>')
      for cell in self.cell:
        cell.write(file,type)
      file.write('</tr>\n')

# You shouldn't need to use this class - best to use Table
class Cell(Content):
  """
  Class to store a table cell
  """
  def __init__(self,root):
    self.contentList = []
    self.table = []
    self.root = root

  def write(self, file, type):
    if type == 'IUL':
      file.write('<td>')
      for content in self.contentList:
        content.write(file,type)
      file.write('</td>')


# finally something that does something. (but not very much)
class Link:

  def __init__(self, link, text):
    self.link = link
    self.text = text

  def write(self, file, type, number=0):
    if type == 'IUL':
      file.write('<a href="' + self.link + '">' + self.text + '</a>')

    if type == 'talkBack':
      file.write('['+str(number)+'.Link]\n')
      file.write('link='+self.link+'\n')
      file.write('text='+self.text+'\n')
      
class Verbatim:

  def __init__(self, text):
    self.text = text

  def write(self, file, type,number=0):
    if type == 'IUL':
      file.write('\n<br><pre>' + self.text + '</pre><br>\n')
    if type == 'talkBack':
      file.write('['+str(number)+'.Verbatim]\n')
      file.write('text='+self.text+'\n')

class Text:

  def __init__(self,text,type='',color=''):
    self.text = text
    self.color = color
    self.type = type

  def write(self, file, type,number =0):
    if type == 'IUL':
      file.write('<p>')
      if self.color:
        file.write('<font color=' + self.color + '>')
      if self.type:
        if isinstance(self.type, str):
          file.write('<'+self.type+'>')
        if isinstance(self.type, list):
          for i in self.type:
            file.write('<'+i+'>')
      file.write(self.text)
      if self.type:
        if isinstance(self.type, str):
          file.write('</'+self.type+'>')
        if isinstance(self.type, list):
          for i in self.type:
            file.write('</'+i+'>')
      if self.color:
        file.write('</font>')
      file.write('</p>\n')

    if type==('talkBack'): 
      file.write('['+str(number)+'.Text]\n')
      if self.type:
        if isinstance(self.type, str):
          file.write('type='+self.type+'\n')
        if isinstance(self.type, list):
          iStr = ''
          for i in self.type:
            iStr += i+',' 
          file.write('type='+iStr[0:-1]+'\n')
      if self.color: file.write('color='+self.color+'\n') 
      file.write('text='+self.text +'\n')

class List:
 
  def __init__(self,list):
    # valid objects in this list are lists of Link, Verbatim, Text, List, Image
    # or any custom object with a compliant write method.
    # remember this expects a list of lists!
    self.list = list
  
  def write(self, file, type,number=0):
    if type == 'IUL':
      file.write('<ol>\n')
      for item in self.list:
        file.write('<li>')
        for element in item:
          element.write(file,type)
        file.write('\n')
      file.write('</ol>\n')
      
class Image:

  def __init__(self,image,root,link=None,width=400):
    self.link = link
    self.image = root+image
    self.width = width
  def write(self, file, type, number = 0):
    if type == 'IUL':
      if self.link:
        file.write('<a href="'+self.link+'"><img src="'+self.image+'" width = '+str(self.width)+'></a>')
      else: 
        file.write('<img src="'+self.image+'" width = '+str(self.width)+'>')   

    if type == 'talkBack':
      file.write('['+str(number)+'.Image]\n')
      file.write('image='+self.image+'\n')
      if self.link: file.write('link='+self.link+'\n')
      if self.width: file.write('width='+self.width+'\n')

    
class Break:

  def __init__(self, times = 1):
    self.times = times
    self.timesList = range(times)

  def write(self, file, type,number=0):
    if type == 'IUL':
      for time in self.timesList:
        file.write('<br>')
    if type == 'talkBack':
      file.write('['+str(number)+'.Linebreak]\n')
      file.write('times='+str(self.times)+'\n')


    
