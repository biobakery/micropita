Installation instructions for microPITA  in a galaxy environment

1. Create a micropita directory under "galaxy-dist/tools"
2. Populate that directory with the following members:
	micropita_input.xml
	micropita.xml
	micropita_prepare.py  
	upload.py
	micropita_format_input_selector.py  
	MicroPITA.py
3. Update member tool_conf.xml  in the galaxy directory adding the following: 
  <section name="micropita" id="micropita">
    <tool file="micropita/micropita_input.xml"/>
    <tool file="micropita/micropita.xml"/>
  </section>
4. Update member datatypes_conf.xml  in the galaxy directory adding the following:
	<datatype extension="micropita" type="galaxy.datatypes.data:Text" subclass="true" display_in_upload="true"/>
5. Recycle galaxy

