Installation instructions for microPITA in a galaxy environment.
These instructions require the Mercurial versioning system, galaxy, and an internet connection.

1. In the  "galaxy-dist/tools" directory install micropita by typing in a terminal:
hg clone https://bitbucket.org/timothyltickle/micropita

2. Update member tool_conf.xml  in the galaxy directory adding the following: 
  <section name="micropita" id="micropita">
    <tool file="micropita/galaxy/micropita_input.xml"/>
    <tool file="micropita/galaxy/micropita.xml"/>
  </section>

3. Update member datatypes_conf.xml  in the galaxy directory adding the following:
	<datatype extension="micropita" type="galaxy.datatypes.data:Text" subclass="true" display_in_upload="true"/>

4. Copy member HMPStool10PCoA.png  to /galaxy/static/images/HMPStool10PCoA.png

5. Recycle galaxy

