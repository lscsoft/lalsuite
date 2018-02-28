<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns="http://www.w3.org/1999/xhtml" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:template match="LIGO_LW">
	<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
		<head>
		    <style type="text/css">
				body {
				    font-family: Helvetica, sans-serif;
				    font-size: 11px;
					text-align: left;
				}

				h2 {
				    font-size: 13px;
				    margin: 10px 0px 3px 0px;
				}

				.parametertable {
					width: 1024px;
				}

				table {
				    width: 1500px;
				    border: 1px solid #444;
				    border-collapse: collapse;
				}

				tr {
					border: 1px dotted #444;
				}
				
				td {
					border: 1px dotted #444;
				}
			</style>
		    <title>LIGO_LW File Format</title>
		</head>
		<body>
			<xsl:for-each select="Table">
				<h2><xsl:value-of select="@Name"/></h2>
				<div class="parametertable"><table>
					<tr><xsl:call-template name="handle-table-columns"/></tr>
					<xsl:for-each select="Stream">
						<xsl:call-template name="handle-table-lines"/>
					</xsl:for-each>
				</table>
				<!-- Check if there's a figure... -->
            	<xsl:if test="contains(./@Name,'summ_mimegroup:summ_mime:table')">
                    <xsl:if test="contains(./Stream/.,'png')">
                        <xsl:call-template name="get-fig-name">
                            <xsl:with-param name="figname" select="normalize-space(substring-before(./Stream/.,'.png'))"/>
                        </xsl:call-template>
                    </xsl:if>
            	</xsl:if>				
				</div>
			</xsl:for-each>
		</body>
	</html>
</xsl:template>

<xsl:template name="get-fig-name">
    <xsl:param name="figname" select="." />
	<xsl:choose>
		<xsl:when test="contains($figname, '&quot;')">
			<xsl:call-template name="get-fig-name">
				<xsl:with-param name="figname" select="substring-after($figname, '&quot;')"/>
			</xsl:call-template>
		</xsl:when>
		<xsl:otherwise>
			<img src="{$figname}.png"/>
		</xsl:otherwise>
	</xsl:choose>    
</xsl:template>

<xsl:template name="handle-table-columns">
	<xsl:for-each select="Column">
		<td><b>
			<xsl:value-of select="position()"/>: 
			<xsl:call-template name="drop-namespace">
				<xsl:with-param name="string" select="./@Name"/>
			</xsl:call-template>
		</b></td>
	</xsl:for-each>
</xsl:template>

<xsl:template name="drop-namespace">
	<xsl:param name="string" select="."/>
	<xsl:choose>
		<xsl:when test="contains($string,':')">
			<xsl:call-template name="drop-namespace">
				<xsl:with-param name="string" select="substring-after($string,':')"/>
			</xsl:call-template>
		</xsl:when>
		<xsl:otherwise>
			<xsl:value-of select="$string"/>
		</xsl:otherwise>
	</xsl:choose>
</xsl:template>

<xsl:template name="handle-table-lines">
	<xsl:param name="string" select="." />
	<xsl:choose>
		<xsl:when test="contains($string, '&#xA;')">
			<tr>
				<xsl:call-template name="divide">
					<xsl:with-param name="to-be-divided" select="normalize-space(substring-before($string, '&#xA;'))" />
				</xsl:call-template>
			</tr>
			<xsl:call-template name="handle-table-lines">
				<xsl:with-param name="string" select="substring-after($string, '&#xA;')" />
			</xsl:call-template>
		</xsl:when>
		<xsl:otherwise>
			<tr>
				<xsl:call-template name="divide">
					<xsl:with-param name="to-be-divided" select="normalize-space($string)" />
				</xsl:call-template>
			</tr>
		</xsl:otherwise>
	</xsl:choose>
</xsl:template>

<xsl:template name="divide">
<xsl:param name="to-be-divided"/>
<xsl:if test="$to-be-divided">
	<xsl:choose>
		<xsl:when test="starts-with($to-be-divided,'&#x22;')">
			<!-- if the string start with '"' -->
			<xsl:choose>
				<xsl:when test="contains(substring-after($to-be-divided,'&#x22;'),'&#x22;,')">
					<!-- if the string starts with '"', see if there is a matching '",' to end it, and recurse -->
					<td><xsl:value-of select="substring-before(substring-after($to-be-divided,'&#x22;'),'&#x22;,')"/></td>
					<xsl:call-template name="divide">
				    	<xsl:with-param name="to-be-divided" select="substring-after(substring-after($to-be-divided,'&#x22;'),'&#x22;,')"/>
					</xsl:call-template>
				</xsl:when>
				<xsl:otherwise>
					<!-- otherwise end at the matching '"', and stop recursion -->
					<td><xsl:value-of select="substring-before(substring-after($to-be-divided,'&#x22;'),'&#x22;')"/></td>										
				</xsl:otherwise>
			</xsl:choose>
		</xsl:when>
		<xsl:otherwise>
			<!-- if the string does not start with '"' -->
			<xsl:choose>
				<xsl:when test="contains($to-be-divided,',')">
					<!-- if there is a comma, split before it, and recurse -->
					<td><xsl:value-of select="substring-before($to-be-divided,',')"/></td>
					<xsl:call-template name="divide">
				    	<xsl:with-param name="to-be-divided" select="substring-after($to-be-divided,',')"/>
					</xsl:call-template>
				</xsl:when>
				<xsl:otherwise>
					<!-- otherwise stop recursion -->
					<td><xsl:value-of select="$to-be-divided"/></td>										
				</xsl:otherwise>
			</xsl:choose>
		</xsl:otherwise>	
	</xsl:choose>
</xsl:if>
</xsl:template>
                        
</xsl:stylesheet>
