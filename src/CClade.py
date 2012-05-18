
import blist
import sys

class CClade:
	
        #Initialize CCClade
        #Dictionary to hold children nodes in tree
        #adValues is a list of the abundance values from the samples
	def __init__( self ):
		
		self.m_hashChildren = {}
		self.m_adValues = None
	
        #Recursively travel the length of a tree until you find the terminal node
        #(where astrClade == Falseor actually [])
        #or a dict key that matches the clade call.
        #If at any time a clade is given that is not currently know, return a new clade
        #which is set to the current Clade as a child.
	def get( self, astrClade ):
		
		return self.m_hashChildren.setdefault(
			astrClade[0], CClade( ) ).get( astrClade[1:] ) if astrClade else self

        #Set all the values given as a list in the same order given
	def set( self, adValues ):
		
		self.m_adValues = blist.blist( [0] ) * len( adValues )
		for i, d in enumerate( adValues ):
			if d:
				self.m_adValues[i] = d

        #This allows you to recursively impute values for clades without values given their children counts.
        #Assumably this should be called only once and after all clade abundances have been added.
        #If the m_adValues already exist return the stored m_adValues. (No imputation needed).
        #Otherwise call impute for all children and take the sum of the values from all the children by column
        #Not a sum of a list but summing a list with lists by element.
	def impute( self ):
		
                #If values do not exist
		if not self.m_adValues:
                        #Call impute on all children
                        #If the parent clade has no abundance values
                        #Then take a copy of the child's
                        #If they now have a copy of a child's but have other children
                        #Sum their children with thier current values
			for pChild in self.m_hashChildren.values( ):
				adChild = pChild.impute( )
				if self.m_adValues:
					for i in range( len( adChild or [] ) ):
						if adChild[i]:
							self.m_adValues[i] += adChild[i]
				elif adChild:
					self.m_adValues = adChild[:] 
		#If values exist return			
		return self.m_adValues
	
        #Update the given hashValues dict with clade abundances given depth specifications
        #Return a set of integers returning an indicator of the structure of the tree preserved in the dict/hash
        #When the appropriate level of the tree is reached
        #Hashvalue is updated with the clade (including lineage) and the abundance. looks like {"clade":blist(["0.0","0.1"...])}
	def _freeze( self, hashValues, iTarget, astrClade, iDepth, fLeaves ):
		
                #fHit is true on atleast one of the following conditions:
                #iTarget is not 0 indicating no changes
                #Leaves are indicated to be only given and the target depth for the leaves is reached.
                #The required depth is reached.
		fHit = ( not iTarget ) or ( ( fLeaves and ( iDepth == iTarget ) ) or ( ( not fLeaves ) and ( iDepth <= iTarget ) ) )
                #Increment depth
		iDepth += 1
                #Returns a set
		setiRet = set()
                #If there are children build integer set indicating structure of the tree preserved in the dict
		if self.m_hashChildren:
                        #Union all the results from freeze of all children
                        #Call freeze but append the child clade to the clade in the call.
                        #And give an incremented depth
			for strChild, pChild in self.m_hashChildren.items( ):
				setiRet |= pChild._freeze( hashValues, iTarget, astrClade + [strChild], iDepth, fLeaves )
			setiRet = set( ( i + 1 ) for i in setiRet )
		else:
			setiRet.add( 0 )
                #Indicate if the correct level is reached
		if iTarget < 0:
			if fLeaves:
				fHit = -( iTarget + 1 ) in setiRet
			else:
				fHit = -( iTarget + 1 ) <= max( setiRet )
                #if astrClade is not == [] (so you are actually in a clade of the tree)
                #And the clade has values (should be true, if not impute should have been callded before running this method)
                #And we are at the correct level of the tree then
                #Add to the dict the clade and the abundance values
		if astrClade and self.m_adValues and fHit:
			hashValues["|".join( astrClade )] = self.m_adValues
		return setiRet
	
        #Call helper function setting the clade and depth to defaults (start positions)
        #The important result of this method is hashValues is updated with clade and abundance information
	def freeze( self, hashValues, iTarget, fLeaves ):
		
		self._freeze( hashValues, iTarget, [], 0, fLeaves )
	
        #Represent tree clade for debugging
	def _repr( self, strClade ):

		strRet = "<"
		if strClade:
			strRet += "%s %s" % (strClade, self.m_adValues)
			if self.m_hashChildren:
				strRet += " "
		if self.m_hashChildren:
			strRet += " ".join( p._repr( s ) for (s, p) in self.m_hashChildren.items( ) )
		
		return ( strRet + ">" )
		
	def __repr__( self ):
		
		return self._repr( "" )

"""
pTree = CClade( )
pTree.get( ("A", "B") ).set( [1, 2, 3] )
pTree.get( ("A", "C") ).set( [4, 5, 6] )
pTree.get( ("D", "E") ).set( [7, 8, 9] )
iTaxa = 0
if iTaxa:
	pTree.impute( )
hashFeatures = {}
pTree.freeze( hashFeatures, iTaxa )
print( pTree )
print( hashFeatures )
sys.exit( 0 )
#"""
