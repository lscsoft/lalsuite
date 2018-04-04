<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns="http://www.w3.org/1999/xhtml" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:output method="html" indent="yes" doctype-system="http://www.w3.org/TR/html4/strict.dtd" doctype-public="-//W3C//DTD HTML 4.01//EN" />

<xsl:template match="LIGO_LW">
	<html>
		<head>
			<style type="text/css">
				body { text-align: left; font-family: Helvetica, sans-serif; font-size: 10pt; }
				
				div { padding: 5pt; }
                div.toc { border: 1px solid #f00; }
                div.table { border: 1px solid #000; margin-top: 10pt; }
								
				h2 { margin-top: 0pt; margin-bottom: 5pt; }
				li { font-size: 11pt; }
				
				table { border: 0px; border-collapse: collapse; border-spacing: 0px; }
                th { border: 1px #69C solid; background-color: #aaa; color: #000; padding: 4px; }
                tr { background-color: #fff; }
                tr.blue { background-color: #eef; }             
                td { border: 1px #69C solid; color: #333; padding: 4px; }

				img { border: 0px; }
			</style>
			
			<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.3.2/jquery.min.js"></script>
			<script type="text/javascript" src="ligolw.js"></script>
			<script>
			    $(document).ready(function() { setup(); });
		    </script>

		    <title>LIGO_LW CBC Document</title>
		</head>
		
		<body>
		    <!-- Build the TOC from the names of the tables -->
		    <div class="toc" id="mytoc">
		        <h2>Table of contents</h2>
                <ul><xsl:apply-templates mode="toc"/></ul>
			</div>

            <!-- Render each table -->
            <xsl:apply-templates mode="render"/>
        </body>
    </html>
</xsl:template>

<!-- Table template for toc -->
<xsl:template match="Table" mode="toc">
    <li><xsl:value-of select="@Name"/></li>
</xsl:template>

<!-- Table template for body -->
<xsl:template match="Table" mode="render">
    <div class="table">
        <a name="{@Name}"><h2 class="title"><xsl:value-of select="@Name"/></h2></a>
        
        <table class="content">
            <!-- Handle the column names -->
            <tr>
                <xsl:apply-templates select="Column"/>
            </tr>
        </table>

		<!-- render the stream as a preformatted block -->
		<pre>
			<xsl:value-of select="Stream"/>
		</pre>
    </div>
</xsl:template>

<!-- Column within Table -->
<xsl:template match="Column">
    <th><xsl:value-of select="@Name"/> (<xsl:value-of select="@Type"/>)</th>
</xsl:template>

</xsl:stylesheet>
