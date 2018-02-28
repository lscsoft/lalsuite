<?xml version="1.0"?>

<xsl:stylesheet version="1.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:template match="/">
  <html>
      <body>
      <xsl:for-each select="monitor">
       
      <h1><xsl:value-of select="name"/> Monitor! </h1>
      <h2><xsl:value-of select="name"/> Events</h2>
        <table border="1">
            <tr bgcolor="#9acd32">
                <th>Type</th>
                <th>Ifo</th>
                <th>Name</th>
                <th>Description</th>
                <th>gps_time</th>
                <th>latency</th>
                <th>Good/Bad Boolean</th>
                <th>Process Duration GPS Difference</th>
                <th>Process Duration Boolean</th>
            </tr>
            <xsl:for-each select="event">
               <tr>
                 <td><xsl:value-of select="type"/></td>
                 <td><xsl:value-of select="ifo"/></td>
                 <td><xsl:value-of select="name"/></td>
                 <td><xsl:value-of select="description"/></td>
                 <td><xsl:value-of select="gps_time"/></td>
                 <td><xsl:value-of select="latency"/></td>
                 <td><xsl:value-of select="bool"/></td>
                 <td><xsl:value-of select="process_duration/gps_difference"/></td>
                 <td><xsl:value-of select="process_duration/bool"/></td>
               </tr>
             </xsl:for-each>
        </table>
      <h2><xsl:value-of select="name"/> Process Durations</h2>
        <table border="1">
            <tr bgcolor="#9acd32">
                <th>IFO</th>
                <th>Name</th>
                <th>Description</th>
                <th>Good/Bad Boolean</th>
                <th>Duration GPS Difference</th>
                <th>Boolean</th>
            </tr>
            <xsl:for-each select="process_duration">
               <tr>
                 <td><xsl:value-of select="ifo"/></td>
                 <td><xsl:value-of select="name"/></td>
                 <td><xsl:value-of select="description"/></td>
                 <td><xsl:value-of select="latency"/></td>
                 <td><xsl:value-of select="gps_difference"/></td>
                 <td><xsl:value-of select="bool"/></td>
               </tr>
             </xsl:for-each>
        </table>
      <h2><xsl:value-of select="name"/> Latency Comparisons</h2>
        <table border="1">
            <tr bgcolor="#9acd32">
                <th>Event 1 Name</th>
                <th>Event 1 Description</th>
                <th>Event 2 Name</th>
                <th>Event 2 Description</th>
                <th>Acceptable Window</th>
                <th>GPS Difference</th>
            </tr>
            <xsl:for-each select="latency_comparison">
               <tr>
                   <xsl:for-each select="event">
                      <td><xsl:value-of select="event/name"/></td>
                      <td><xsl:value-of select="event/description"/></td>
                   </xsl:for-each>
                   <td><xsl:value-of select="window"/></td>
                   <td><xsl:value-of select="gps_difference"/></td>
               </tr>
             </xsl:for-each>
         </table>
     </xsl:for-each>
  </body>
  </html>
</xsl:template>
</xsl:stylesheet> 

