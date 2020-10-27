/** \file xmlfile.h
  * \brief XML input file
  * \author Martin Bernreuther <bernreuther@hlrs.de>

  * encapsulates XML main input file
  * handles include functionality
  * query syntax similar to (subset of) XPATH

  * based on rapidxml (http://rapidxml.sourceforge.net/)
  * similar functionality can be found in the Boost Property tree class
    (http://www.boost.org/doc/libs/1_41_0/doc/html/property_tree.html)
  * use XPATH, if XML library supports it.
    (like e.g. tinyxpath(http://tinyxpath.sourceforge.net/), which is a successor of
               tinyxml(http://www.grinninglizard.com/tinyxml/) and might one day also included into
               TinyXML++ (http://code.google.com/p/ticpp/)
     or a "heavier" library like the Xerces-C++ (http://xerces.apache.org/xerces-c/) based
               Xalan-C++(http://xml.apache.org/xalan-c/) or XQilla(http://xqilla.sourceforge.net/)
            or libxml2 (http://xmlsoft.org/)
    )
*/

#ifndef XMLFILE_H
#define XMLFILE_H

#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <set>
#include <stack>

//#define RAPIDXML_NO_EXCEPTIONS
#include "rapidxml/rapidxml.hpp"

#ifdef ENABLE_MPI
#include "mpi.h"
#endif

//+XMLfile======================================================================
/**
* \class XMLfile
* \brief XML file abstraction
*
* DOM representation of an XML file
*/
class XMLfile
{
	//typedef rapidxml::xml_node<> t_XMLnode;
	typedef rapidxml::xml_base<> t_XMLnode;
	typedef rapidxml::xml_node<> t_XMLelement;
	typedef rapidxml::xml_attribute<> t_XMLattribute;
	typedef rapidxml::xml_document<> t_XMLdocument;
	
public:
	static const char *const includeattrtag;
	static const char *const queryattrtag;

	class Query;
	//class Query::const_iterator; forward declaration not possible...
	
	
//+XMLfile::Node----------------------------------------------------------------
	class Node
	{
	friend class Query;
	friend class XMLfile;
	// using member function friends is more specific about where private access happens, but is difficult to handle for a one-pass-compiler...
	//friend Node Query::operator [](unsigned long idx) const;
	//friend Node Query::const_iterator::operator *() const;
	//friend bool changecurrentnode(const Query::const_iterator& pos);
	//friend unsigned long changecurrentnode(const std::string& nodepath);
	//friend void expandincludes();
	//friend unsigned long query(std::list<Node>& nodeselection, const char* querystr, Node startnode) const;
	//friend bool initfile_local(const std::string& filepath);
	//friend void initstring_local(const char* xmlstring);
	
	public:
		/// \enum NodeType
		/// \brief Node types enumeration
		/// associates a number to each XML node type
		/// \note cmp. http://www.w3schools.com/Dom/dom_nodetype.asp
		enum NodeType { Unknown_Node=-1,
		                 Invalid_Node=0,
		                 ELEMENT_Node=1,
		                 ATTRIBUTE_Node=2,
		                 TEXT_Node=3,
		                 CDATA_SECTION_Node=4,
		                 ENTITY_REFERENCE_Node=5,
		                 ENTITY_Node=6,
		                 PROCESSING_INSTRUCTION_Node=7,
		                 COMMENT_Node=8,
		                 DOCUMENT_Node=9,
		                 DOCUMENT_TYPE_Node=10,
		                 DOCUMENT_FRAGMENT_Node=11,
		                 NOTATION_Node=12 };
		
		/// \brief XMLfile::Node constructor
		/// sets up an invalid node
		Node()
			: m_xmlnode(NULL), m_nodepath(std::string()), m_type(Unknown_Node)
			{}
		
		/// \brief XMLfile::Node copy constructor
		/// duplicate a node
		/// \param	XMLfile::Node&	source node
		Node(const Node& n)
			: m_xmlnode(n.m_xmlnode), m_nodepath(n.m_nodepath), m_type(n.m_type)
			{}
		
		/// \brief get nodepath
		/// returns the full path of the node
		/// \return	std::string	nodepath
		const std::string& nodepath() const
			{ return m_nodepath; }
		
		/// \brief get node name
		/// returns the name of the node, if node is valid - otherwise an empty string
		std::string name() const
			{ return (m_xmlnode) ? std::string(m_xmlnode->name()) : std::string(); }
		
		/// \brief invalidate the node
		/// set the internal state to invalid
		void invalidate()
			{ m_xmlnode=NULL; m_nodepath=std::string(); }
		
		/// \brief get node type
		/// determine the node type (enumeration)
		/// \return	XMLfile::Node::NodeType	node type
		NodeType type() const
			{ return m_type; }
		
		/// \brief check, if node is root node 
		/// return true, if node is a valid root node - false otherwise
		/// \return	bool
		bool isRootNode() const
			{ if(m_xmlnode && ! m_xmlnode->parent()) return true; else return false; }
		/// \brief check, if node is leaf node 
		/// return true, if node is a valid leaf node - false otherwise
		/// \return	bool
		bool isLeafNode() const;
		
		/// \brief get the node value
		/// returns the node value converted to a given type
		/// \param	T&	Value to return
		/// \return	bool success?
		template<typename T> bool getValue(T& value) const;
		
		/// \brief get the node string value
		/// returns the node value as string value or default value, if node is invalid
		/// \param	std::string	default value (default: empty string)
		/// \return	std:string	return value
		std::string value_string(std::string defaultvalue=std::string()) const;
		/// \brief get the node int value
		/// returns the node value as int value or default value, if node is invalid
		/// \param	int	default value (default: 0)
		/// \return	int	return value
		int value_int(int defaultvalue=0) const;
		/// \brief get the node long value
		/// returns the node value as long value or default value, if node is invalid
		/// \param	long	default value (default: 0)
		/// \return	long	return value
		long value_long(long defaultvalue=0) const;
		/// \brief get the node float value
		/// returns the node value as float value or default value, if node is invalid
		/// \param	float	default value (default: 0.)
		/// \return	float	return value
		float value_float(float defaultvalue=0.) const;
		/// \brief get the node double value
		/// returns the node value as double value or default value, if node is invalid
		/// \param	double	default value (default: 0.)
		/// \return	double	return value
		double value_double(double defaultvalue=0.) const;
		/// \brief get the node bool value
		/// returns the node value as bool value or default value, if node is invalid
		/// \param	bool	default value (default: false)
		/// \return	bool	return value
		bool value_bool(bool defaultvalue=false) const;
		
		/// \brief assignment operator
		/// copy/duplicate other node content to node
		/// \param	XMLfile::Node	source node
		/// \return	XMLfile::Node&	reference to this node
		Node& operator =(const Node& n)
			{ m_xmlnode=n.m_xmlnode; m_nodepath=n.m_nodepath; return *this; }
		/// \brief template type cast operator 
		/// return node content (getValue)
		template<typename T> operator T() const
			{ T value; getValue(value); return value; }
		/// \brief bool cast operator 
		/// check if Node is valid
		operator bool() const
			{ return m_xmlnode!=NULL && m_type!=Invalid_Node; }
		
		/// \brief print XML data to stream
		/// print the node content using XML syntax
		/// \param	std::ostream&	stream to write to (default: std::cout)
		void printXML(std::ostream& ostrm=std::cout) const;
		/// \brief print data to stream
		/// print the node data
		/// \param	std::ostream&	stream to write to (default: std::cout)
		void print(std::ostream& ostrm=std::cout) const;
	
		Node& operator=(Node& node){
			m_xmlnode = node.m_xmlnode;
			m_nodepath = node.m_nodepath;
			m_type = node.m_type;
			return *this;
		}
	private:
		Node(const t_XMLnode* xmlnode, std::string nodepath=std::string())
		    : m_xmlnode(xmlnode), m_nodepath(nodepath), m_type(Unknown_Node)
		    { if(!xmlnode) m_type=Invalid_Node; }
		Node(const t_XMLelement* xmlelement, std::string nodepath=std::string());
		Node(const t_XMLattribute* xmlattribute, std::string nodepath=std::string())
		    : m_xmlnode(xmlattribute), m_nodepath(nodepath), m_type(ATTRIBUTE_Node)
		    { if(!xmlattribute) m_type=Invalid_Node; }
		/* rapidxml::xml_document is derived from rapidxml::xml_node => t_XMLelement* constructor will handle this
		Node(const t_XMLdocument* xmldocument, std::string nodepath=std::string())
		    : m_xmlnode(xmldocument), m_nodepath(nodepath), m_type(DOCUMENT_Node)
		    { if(!xmldocument) m_type=Invalid_Node; }
		*/
		
		const t_XMLnode*	m_xmlnode;
		std::string	m_nodepath;
		NodeType	m_type;
	};
//-XMLfile::Node----------------------------------------------------------------
	
	
//+XMLfile::Query---------------------------------------------------------------
	class Query
	{
	friend class XMLfile;
	//friend Query XMLfile::query(const char* querystr) const;
	
	public:
	
//+XMLfile::Query::const_iterator...............................................
		class const_iterator
		{
		friend class Query;
	
		public:
			/// \brief XMLfile::Query::const_iterator constructor
			/// sets up an invalid iterator
			const_iterator()
				: m_query(NULL), m_nodesidx(0) {}
			long index() const { return m_nodesidx; }
			/// \brief bool operator
			/// check if iterator is valid and index is in a valid range
			operator bool() const
				{ return m_query!=NULL && m_nodesidx>=0 && m_nodesidx<long(m_query->card()); }
			/// \brief prefix increment operator
			/// put the index to the next position
			/// \return XMLfile::Query::const_iterator&	reference to this iterator
			const_iterator& operator ++()
				{ if(m_nodesidx<long(m_query->card())) ++m_nodesidx; return *this; }
			/// \brief postfix increment operator
			/// increment the index
			/// \return XMLfile::Query::const_iterator&	reference to this iterator
			const_iterator operator ++(int)
				{ if(m_nodesidx<long(m_query->card())) m_nodesidx++; return *this; }
			/// \brief prefix decrement operator
			/// decrement the index
			/// \return XMLfile::Query::const_iterator&	reference to this iterator
			const_iterator& operator --() // prefix decrement
				{ if(m_nodesidx>=0) --m_nodesidx; return *this; }
			/// \brief postfix decrement operator
			/// decrement the index
			/// \return XMLfile::Query::const_iterator&	reference to this iterator
			const_iterator operator --(int)     // postfix decrement
				{ if(m_nodesidx>=0) m_nodesidx--; return *this; }
			/// \brief access operator
			/// return the node, the iterator (index) actually points to or an invalid Node, if the iterator is invalid
			/// \return XMLfile::Node	actual Node
			Node operator *() const
				{ if(*this) return (*m_query)[m_nodesidx]; else return Node(); }
	
		private:
			const Query* m_query;
			long m_nodesidx;
	
			const_iterator(const Query* query, long idx)
				: m_query(query), m_nodesidx(idx) {}
		};
//-XMLfile::Query::const_iterator...............................................
		
		/// \brief XMLfile::Query constructor
		/// sets up an invalid query
		Query() : m_xmlfile(NULL), m_nodes()
			{}
		/// \brief copy constructor
		/// duplicate a given query and register the new created one at the XMLfile
		/// \param const Query&	query
		Query(const Query& q) : m_xmlfile(q.m_xmlfile), m_nodes(q.m_nodes)
			{ xmlfile_register(); }
		/// \brief XMLfile::Query destructor
		/// unregister the query
		~Query()
			{ xmlfile_unregister(); }
		
		/// \brief get cardinality of query set
		/// return the size of the query set
		/// \return unsigned long	cardinality
		unsigned long card() const { return m_nodes.size(); }
		/// \brief check ,if query set is empty
		/// return true if the query set is empty
		/// \return bool	is empty?
		bool empty() const
			{ return card()==0; }
		/// \brief invalidate query
		/// unregister from XMLfile and invalidate the query
		void invalidate()
			{ xmlfile_unregister(); m_xmlfile=NULL; m_nodes.clear(); }
		/// \brief get first query entry
		/// return the first Node of the query set (invalid Node, if there's none)
		/// \return XMLfile::Node	first Node
		Node front() const
			{ return (*this)[0]; }
		
		/// \brief assignment operator
		/// copy/duplicate other query content to query
		/// \param const Query&	source query
		/// \return XMLfile::Query	reference to this query
		Query& operator =(const Query& q);
		/// \brief unsigned long cast operator
		/// return cardinality
		operator unsigned long() const
			{ return card(); }
		/// \brief bool cast operator
		/// return if query is valid
		operator bool() const
			{ return m_xmlfile!=NULL; }
		/// \brief indexing operator
		/// return the Node within the query at given index, or empty node, if query or index is invalid
		/// \param unsigned long	index
		/// \return XMLfile::Node	node
		Node operator [](unsigned long idx) const
			{ if(idx<m_nodes.size()) return Node(m_nodes[idx]); else return Node(); }
		
		/// \brief get starting iterator
		/// return an iterator to the first node
		/// \return XMLfile::Query::const_iterator	iterator
		const_iterator begin() const
			{ return const_iterator(this,0); }
		/// \brief get reverse starting iterator
		/// return an iterator to the last node
		/// \return XMLfile::Query::const_iterator	iterator
		const_iterator rbegin() const
			{ return const_iterator(this,m_nodes.size()-1); }
		/// \brief get (dummy) iteration end for comparison
		/// return false (as dummy) to enable e.g. for(it=query.begin();it!=query.end();++it)
		/// this is preferable to just use e.g. for(it=query.begin();it;++it)
		/// \return boolean	false
		bool end() const
			{ return false; }
		/*
		/// \brief get (dummy) iteration end for comparison
		/// return false (as dummy) for enable e.g. for(it=query.rbegin();it!=query.rend();--it)
		/// Note: A reverse iterator should proceed with ++it (and therefore probably a "reverse" flag iterator attribute is needed)
		/// \return boolean	false
		bool rend() const
			{ return false; }
		*/
		
		/// \brief get node value
		/// get the node content and convert it to a given type
		/// \param T&	variable to return value
		/// \return unsigned long	number of nodes matching the nodepath
		template<typename T> unsigned long getNodeValue(T& value) const;
		/// \brief get node value as string
		/// get the node content
		/// \param const std::string	default value to return, if node is not found
		/// \return std::string	node value
		std::string getNodeValue_string(const std::string defaultvalue=std::string()) const
			{ std::string value(defaultvalue); getNodeValue(value); return value; }
		/// \brief get node value as int
		/// get the node content and convert it to an int
		/// \param int	default value to return, if node is not found
		/// \return int	node value
		int getNodeValue_int(int defaultvalue=0) const
			{ int value=defaultvalue; getNodeValue(value); return value; }
		/// \brief get node value as long
		/// get the node content and convert it to a long
		/// \param long	default value to return, if node is not found
		/// \return long	node value
		long getNodeValue_long(long defaultvalue=0) const
			{ long value=defaultvalue; getNodeValue(value); return value; }
		/// \brief get node value as float
		/// get the node content and convert it to a float
		/// \param float	default value to return, if node is not found
		/// \return float	node value
		float getNodeValue_float(float defaultvalue=0.) const
			{ float value=defaultvalue; getNodeValue(value); return value; }
		/// \brief get node value as double
		/// get the node content and convert it to a double
		/// \param double	default value to return, if node is not found
		/// \return double	node value
		double getNodeValue_double(double defaultvalue=0.) const
			{ double value=defaultvalue; getNodeValue(value); return value; }
		/// \brief get node value as bool
		/// get the node content and convert it to a bool
		/// \param bool	default value to return, if node is not found
		/// \return bool	node value
		bool getNodeValue_bool(bool defaultvalue=false) const
			{ bool value=defaultvalue; getNodeValue(value); return value; }
		
		/// \brief print XML data to stream
		/// print the query content using XML syntax
		/// \param	std::ostream&	stream to write to (default: std::cout)
		void printXML(std::ostream& ostrm=std::cout) const;
		/// \brief print data to stream
		/// print the query content
		/// \param	std::ostream&	stream to write to (default: std::cout)
		void print(std::ostream& ostrm=std::cout) const;
		
	private:
		const XMLfile*	m_xmlfile;
		std::vector<Node> m_nodes;    // nodes within the actual query set
		
		Query(const XMLfile* xmlfile)
			: m_xmlfile(xmlfile)
			{ m_nodes.clear(); }
		
		void xmlfile_register()
			{ if(m_xmlfile) m_xmlfile->registerQuery(this); }
		void xmlfile_unregister()
			{ if(m_xmlfile) m_xmlfile->unregisterQuery(this); }
	};
//-XMLfile::Query---------------------------------------------------------------
	
		friend void Query::xmlfile_register();
		friend void Query::xmlfile_unregister();
	
	/// \brief XMLfile default constructor
	XMLfile();
	/// \brief XMLfile default destructor
	virtual ~XMLfile() { clear(); }
	
	/// \brief constructor for XML-file
	/// constructor calls initfile
	/// \param std::string&	XML-file
	XMLfile(const std::string& filepath);
	
	/// \brief initialize with XML-file
	/// instantiating with XML file
	/// \param std::string&	XML-file
	bool initfile(const std::string& filepath);
	
	/// \brief initialize with XML-string
	/// instantiating with XML-string
	/// \param const char*	XML-string
	void initstring(const char* xmlstring);
	
	/// \brief get XML file directory
	/// if instantiated with a XML-file, return the directory
	/// \return std::string	directory
	const std::string getDir() const { return m_filedir; }
	/// \brief get XML filename
	/// if instantiated with a XML-file, return the filename (without directory part of path)
	/// \return std::string	filename
	const std::string getFilename() const { return m_filename; }
	
	/// \brief set current node
	/// set a node, relative queries start with
	/// \param std::string&	node path (default: "/")
	/// \return long	cardinality of the resulting query set; <=0 => no action
	long changecurrentnode(const std::string& nodepath=std::string("/"));
	/// \brief set current node
	/// set a node, relative queries start with
	/// \param Query::const_iterator&	query iterator pointing to a node
	/// \return bool	done?
	bool changecurrentnode(const Query::const_iterator& pos);
	/// \brief get current node path
	/// 
	/// \return std::string	node path
	std::string getcurrentnodepath() const { return m_currentnode.nodepath(); }
	
	/// \brief get node value
	/// get the node content and convert it to a given type
	/// \param std::string&	nodepath
	/// \param T&	variable to return value
	/// \return unsigned long	number of nodes matching the nodepath
	template<typename T> unsigned long getNodeValue(const std::string& nodepath, T& value) const
		//{ Query q=query(nodepath.c_str()); if(!q.empty()) q.front().getValue(value); return q.card(); }
		{ return query(nodepath).getNodeValue(value); }
	/// \brief get node value
	/// get the node content and convert it to a given type
	/// \param const char*	nodepath
	/// \param T&	variable to return value
	/// \return unsigned long	number of nodes matching the nodepath
	template<typename T> unsigned long getNodeValue(const char* nodepath, T& value) const
		{ return getNodeValue(std::string(nodepath),value); }
	/// \brief get node value as string
	/// get the node content
	/// \param std::string&	nodepath
	/// \param std::string&	default value to return, if node is not found
	/// \return std::string	node value
	std::string getNodeValue_string(const std::string& nodepath, const std::string defaultvalue=std::string()) const
		//{ std::string value(defaultvalue); getNodeValue(nodepath,value); return value; }
		{ return query(nodepath).getNodeValue_string(defaultvalue); }
	/// \brief get node value as string
	/// get the node content
	/// \param const char*	nodepath
	/// \param std::string&	default value to return, if node is not found
	/// \return std::string	node value
	std::string getNodeValue_string(const char* nodepath, std::string defaultvalue=std::string()) const
		{ return getNodeValue_string(std::string(nodepath),defaultvalue); }
	/// \brief get node value as int
	/// get the node content and convert it to an integer
	/// \param std::string&	nodepath
	/// \param int&	default value to return, if node is not found
	/// \return int	node value
	int getNodeValue_int(const std::string& nodepath, int defaultvalue=0) const
		//{ int value=defaultvalue; getNodeValue(nodepath,value); return value; }
		{ return query(nodepath).getNodeValue_int(defaultvalue); }
	/// \brief get node value as int
	/// get the node content and convert it to an integer
	/// \param const char*	nodepath
	/// \param int&	default value to return, if node is not found
	/// \return int	node value
	int getNodeValue_int(const char* nodepath, int defaultvalue=0) const
		{ return getNodeValue_int(std::string(nodepath),defaultvalue); }
	/// \brief get node value as long
	/// get the node content and convert it to a long
	/// \param std::string&	nodepath
	/// \param long&	default value to return, if node is not found
	/// \return long	node value
	long getNodeValue_long(const std::string& nodepath, long defaultvalue=0) const
		//{ long value=defaultvalue; getNodeValue(nodepath,value); return value; }
		{ return query(nodepath.c_str()).getNodeValue_long(defaultvalue); }
	/// \brief get node value as long
	/// get the node content and convert it to a long
	/// \param const char*	nodepath
	/// \param long&	default value to return, if node is not found
	/// \return long	node value
	long getNodeValue_long(const char* nodepath, long defaultvalue=0) const
		{ return getNodeValue_long(std::string(nodepath),defaultvalue); }
	/// \brief get node value as float
	/// get the node content and convert it to a float
	/// \param std::string&	nodepath
	/// \param float&	default value to return, if node is not found
	/// \return float	node value
	float getNodeValue_float(const std::string& nodepath, float defaultvalue=0.) const
		//{ float value=defaultvalue; getNodeValue(nodepath,value); return value; }
		{ return query(nodepath).getNodeValue_float(defaultvalue); }
	/// \brief get node value as float
	/// get the node content and convert it to a float
	/// \param const char*	nodepath
	/// \param float&	default value to return, if node is not found
	/// \return float	node value
	float getNodeValue_float(const char* nodepath, float defaultvalue=0.) const
		{ return getNodeValue_float(std::string(nodepath),defaultvalue); }
	/// \brief get node value as double
	/// get the node content and convert it to a double
	/// \param std::string&	nodepath
	/// \param double&	default value to return, if node is not found
	/// \return double	node value
	double getNodeValue_double(const std::string& nodepath, double defaultvalue=0.) const
		//{ double value=defaultvalue; getNodeValue(nodepath,value); return value; }
		{ return query(nodepath).getNodeValue_double(defaultvalue); }
	/// \brief get node value as double
	/// get the node content and convert it to a double
	/// \param const char*	nodepath
	/// \param double&	default value to return, if node is not found
	/// \return double	node value
	double getNodeValue_double(const char* nodepath, double defaultvalue=0.) const
		{ return getNodeValue_double(std::string(nodepath),defaultvalue); }
	/// \brief get node value as bool
	/// get the node content and convert it to a bool
	/// (an alternative of using <tag>true</tag> is to look for an occurance of an empty element like <tag/>)
	/// \param std::string&	nodepath
	/// \param bool&	default value to return, if node is not found
	/// \return bool	node value
	bool getNodeValue_bool(const std::string& nodepath, bool defaultvalue=false) const
		//{ bool value=defaultvalue; getNodeValue(nodepath,value); return value; }
		{ return query(nodepath).getNodeValue_bool(defaultvalue); }
	/// \brief get node value as bool
	/// get the node content and convert it to a bool
	/// \param const char*	nodepath
	/// \param bool&	default value to return, if node is not found
	/// \return bool	node value
	bool getNodeValue_bool(const char* nodepath, bool defaultvalue=false) const
		{ return getNodeValue_bool(std::string(nodepath),defaultvalue); }
	
	/// \brief print node content as XML
	/// print node content to stream using XML
	/// \param std::ostrm&	output stream
	void printXML(std::ostream& ostrm=std::cout) const;
	/// \brief print node content
	/// print node content and debug information to stream
	/// \param std::ostrm&	output stream
	void print(std::ostream& ostrm=std::cout) const;
	/// \brief save
	/// save node content as XML-file
	/// \param std::string&	output file
	void save(std::string filepath=std::string());
	
	/// \brief perform a query
	/// return a query to a given query expression
	/// \param std::string&	query string
	/// \return XMLfile::Query	query
	Query query(const std::string& querystr) const;
	
	/// \brief std::string cast operator
	/// XMLfile will cast to a string with XML content
	operator std::string() const;
	
	/// \brief number of registered queries
	/// return the number of active queries
	/// \return number of active queries
	size_t numqueries() const
		{ return m_queries.size(); }

#ifdef ENABLE_MPI
	void setMPIdefaults(int mpirootrank=0, MPI_Comm mpicomm=MPI_COMM_WORLD);
#endif

//protected:
private:
	std::string m_filedir;
	std::string m_filename;
	t_XMLdocument m_xmldoc;
	Node m_currentnode;
	mutable std::set<Query*> m_queries;
	//std::stack<Node> m_lastnodes;
	
	void clear();
	bool initfile_local(const std::string& filepath);
	void initstring_local(const char* xmlstring);
	void expandincludes();
	unsigned long query(std::list<Node>& nodeselection, const std::string& querystring, Node startnode=Node()) const;
	void insertcloneelement(const t_XMLelement* src, t_XMLelement* dest_after);
	
	void registerQuery(Query* q) const
		{ if(q) m_queries.insert(q); }
	void unregisterQuery(Query* q) const
		{ m_queries.erase(q); }
	void invalidateQueries()
	{ for(std::set<Query*>::iterator pos=m_queries.begin();pos!=m_queries.end();++pos) (*pos)->invalidate(); m_queries.clear(); }

#ifdef ENABLE_MPI
	int m_mpi_rootrank;
	MPI_Comm m_mpi_comm;
	int m_mpi_myrank;
	//bool m_mpi_isroot;

	bool distributeXMLstring();
#endif

};
//-XMLfile======================================================================

/// \brief write a node to a stream
/// write XML node data to an output stream
/// \param std::ostream&	output stream
/// \param XMLfile::Node&	node
/// \return std::ostream&	output stream
inline std::ostream& operator << (std::ostream& ostrm, const XMLfile::Node& xmlnode)
{
	xmlnode.printXML(ostrm);
	return ostrm;
}

/// \brief write a query to a stream
/// write XML query data to an output stream
/// \param std::ostream&	output stream
/// \param XMLfile::Query&	query
/// \return std::ostream&	output stream
inline std::ostream& operator << (std::ostream& ostrm, const XMLfile::Query& xmlquery)
{
	xmlquery.printXML(ostrm);
	return ostrm;
}

/// \brief write a xmlfile to a stream
/// write XML xmlfile data to an output stream
/// \param std::ostream&	output stream
/// \param XMLfile&	xmlfile
/// \return std::ostream&	output stream
inline std::ostream& operator << (std::ostream& ostrm, const XMLfile& xmlfile)
{
	xmlfile.printXML(ostrm);
	return ostrm;
}

/*
// explicit instantiation
template bool XMLfile::Node::getValue<std::string>(std::string& value)const;
template bool XMLfile::Node::getValue<int>(int& value)const;
template bool XMLfile::Node::getValue<long>(long& value)const;
template bool XMLfile::Node::getValue<float>(float& value)const;
template bool XMLfile::Node::getValue<double>(double& value)const;
template bool XMLfile::Node::getValue<bool>(bool& value)const;
*/

#endif
