#Installation instructions for microPITA in a galaxy environment.
These instructions require the Mercurial versioning system, galaxy, and an internet connection.

#For general reference about microPita please refer to:
```
https://bitbucket.org/biobakery/micropita
```




#Installation Instructions
In the  "galaxy-dist/tools" directory install micropita by typing in a terminal:
```
hg clone https://bitbucket.org/biobakery/micropita
```


Update member tool_conf.xml  in the galaxy directory adding the following: 
```
  <section name="micropita" id="micropita">
    <tool file="micropita/galaxy/micropita.xml"/>
  </section>
```

Update member datatypes_conf.xml  in the galaxy directory adding the following:
```
	<datatype extension="micropita" type="galaxy.datatypes.data:Text" subclass="true" display_in_upload="true"/>
```

Copy the 2 *.png  members   to /galaxy/static/images

Recycle galaxy

