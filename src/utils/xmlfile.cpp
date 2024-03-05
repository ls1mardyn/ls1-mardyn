/** \file xmlfile.cpp
  * \brief XML file
  * \author Martin Bernreuther <bernreuther@hlrs.de>
*/

#include "xmlfile.h"

#include <fstream>
#include <sstream>
#include <map>

#include <algorithm>
#include <cstdlib>
#include <string>

#include "utils/mardyn_assert.h"
#include "rapidxml/rapidxml_print.hpp"
#include "String_utils.h"
#include "Simulation.h"

//#include <cstdio>	// fseek(),fread(); should be included after mpi.h
//#ifdef __linux__
//#include <sys/stat.h>   // stat()
//#endif


using namespace rapidxml;

const char *const XMLfile::includeattrtag = "include";
const char *const XMLfile::queryattrtag = "query";


// XMLfile======================================================================
// public methods

XMLfile::XMLfile()
{
	clear();
#ifdef ENABLE_MPI
	setMPIdefaults();
#endif
}

XMLfile::XMLfile(const std::string& filepath)
{
	clear();
#ifdef ENABLE_MPI
	setMPIdefaults();
#endif
	initfile(filepath);
}

bool XMLfile::initfile(const std::string& filepath)
{
	clear();
	bool status;

#ifdef ENABLE_MPI
	// only root node will read XML file in MPI parallel version
	if(m_mpi_myrank==m_mpi_rootrank)
#endif
	status=initfile_local(filepath);

#ifdef ENABLE_MPI
	status=distributeXMLstring();
#endif

	return status;
}

void XMLfile::initstring(const char* xmlstring)
{
	clear();

#ifdef ENABLE_MPI
	// root node XML string will be used in MPI parallel version
	if(m_mpi_myrank==m_mpi_rootrank)
#endif
	initstring_local(xmlstring);

#ifdef ENABLE_MPI
	distributeXMLstring();
#endif
}

long XMLfile::changecurrentnode(const std::string& nodepath)
{
	std::list<Node> nodes;
	unsigned long foundnodes=query(nodes,nodepath.c_str());
	if(foundnodes && nodes.front().type()==Node::ELEMENT_Node)
	{
		//m_lastnodes.push(m_currentnode);
		m_currentnode=nodes.front();
	} else {
		foundnodes=-foundnodes;
	}
	return foundnodes;
}

bool XMLfile::changecurrentnode(const Query::const_iterator& pos)
{
	if(pos && (*pos).type()==Node::ELEMENT_Node)
	{
		//m_lastnodes.push(m_currentnode);
		m_currentnode=*pos;
		return true;
	}
	else
		return false;
}

void XMLfile::save(std::string filepath)
{
#ifdef ENABLE_MPI
	// only root process will save its data for now
	if(m_mpi_myrank!=m_mpi_rootrank) return;
#endif
	if (filepath.empty()) filepath=m_filedir+m_filename;
	std::ofstream of(filepath.c_str());
	printXML(of);
}

XMLfile::Query XMLfile::query(const std::string& querystr) const
{
	std::list<Node> nodes;
	//nodes.clear();
	query(nodes,querystr);

	XMLfile::Query queryresult(this);
	//copy nodes (std::list<Node>) to queryresult.m_nodes (std::vector<Node>)
	queryresult.m_nodes.resize(nodes.size());
/*
	std::vector<Node>::iterator nodespos2=queryresult.m_nodes.begin();
	for(std::list<Node>::const_iterator nodespos=nodes.begin();nodespos!=nodes.end();++nodespos)
	{
		*nodespos2=*nodespos;
		++nodespos2;
	}
*/
	copy(nodes.begin(), nodes.end(), queryresult.m_nodes.begin());
/*
	////queryresult.m_nodes.clear();
	//queryresult.m_nodes.assign(nodes.begin(), nodes.end());
	queryresult.m_nodes.reserve(nodes.size());
	copy(nodes.begin(), nodes.end(), back_inserter(queryresult.m_nodes));
*/
	return queryresult;
}


XMLfile::operator std::string() const
{
	std::stringstream ss;
	ss << m_xmldoc;
	return ss.str();
}

void XMLfile::printXML(std::ostream& ostrm) const
{
	ostrm << m_xmldoc << std::endl;
}

void XMLfile::print(std::ostream& ostrm) const
{
	ostrm << m_filedir << "," << m_filename << " (" << m_currentnode.nodepath() << ")" << std::endl;
	ostrm << m_xmldoc << std::endl;
}


//------------------------------------------------------------------------------
// private methods

void XMLfile::clear()
{
	m_filedir="";
	m_filename="";
	m_xmldoc.clear();
	invalidateQueries();
}

bool XMLfile::initfile_local(const std::string& filepath) {
	clear();
	auto filepathTrimmed = string_utils::trim(filepath);
	std::size_t fnpos=filepathTrimmed.find_last_of('/');
	if (fnpos!=std::string::npos) {
		m_filedir=std::string(filepathTrimmed).substr(0,fnpos+1);
		m_filename=std::string(filepathTrimmed).substr(fnpos+1);
	}
	else {
		m_filedir=std::string();
		m_filename=std::string(filepathTrimmed);
	}

	//version using ifstream
	std::ifstream fstrm(filepathTrimmed.c_str(),std::ifstream::binary|std::ifstream::ate);
	if(!fstrm) {
		std::cerr << "ERROR opening " << filepathTrimmed << std::endl;
		clear();
		Simulation::exit(1);
	}
	std::ifstream::pos_type filesize=fstrm.tellg();
	fstrm.close(); fstrm.clear();
	char* xmlstr = m_xmldoc.allocate_string(nullptr,static_cast<size_t>(filesize)+1);
	xmlstr[filesize]=0;
	//                          std::ios::binary
	fstrm.open(filepathTrimmed.c_str(),std::ifstream::binary);
	//checking if(!fstrm) again should not be necessary
	fstrm.read(xmlstr,filesize);
	filesize-=fstrm.gcount();
	fstrm.close();
//
	mardyn_assert(filesize == 0);

	m_xmldoc.parse<0>(xmlstr);
	expandincludes();
	m_currentnode=Node(&m_xmldoc,"/");

	//cout << "init file " << m_filedir << ";" << m_filename << std::endl;

	return true;
}

void XMLfile::initstring_local(const char* xmlstring)
{
	clear();

	char* xmlstr = m_xmldoc.allocate_string(xmlstring);
	m_xmldoc.parse<0>(xmlstr);
	expandincludes();
	m_currentnode=Node(&m_xmldoc,"/");
}

void XMLfile::expandincludes()
{
	std::list<xml_node<>* > nodelist;
	nodelist.push_back(m_xmldoc.first_node());
	while(!nodelist.empty())
	{
		xml_node<>* node=nodelist.front();
		nodelist.pop_front();
		if(!node) continue;
		if(std::string(node->name())==std::string(includeattrtag))
		{
			std::string inclfile=node->value();
			// relative path?
			if(inclfile.find_first_of("/~")!=0) inclfile=m_filedir+inclfile;
			xml_attribute<>* attr = node->first_attribute(queryattrtag);
			std::string queryval("/");
			if(attr) queryval=attr->value();
			//cout << "XMLfile::expandincludes DEBUG:\tincluding " << inclfile << "\tquery " << queryval << std::endl;
			XMLfile inclxmldoc;
			inclxmldoc.initfile_local(inclfile);
			std::list<Node> inclnodes;
			std::list<Node>::iterator posinclnodes;
			inclxmldoc.query(inclnodes,queryval.c_str());
			for(posinclnodes=inclnodes.begin();posinclnodes!=inclnodes.end();++posinclnodes)
			{
				const xml_node<>* inclnode=static_cast<const xml_node<>*>((*posinclnodes).m_xmlnode);
				//cout << "XMLfile::expandincludes DEBUG:\tinserting " << inclnode->name() << " (" << (*posinclnodes).nodepath() << ")" << std::endl;
				insertcloneelement(inclnode,node);
			}
			xml_node<>* node_parent=node->parent();
			if(node_parent)
				// removing include node itself
				node_parent->remove_node(node);
		}
		for (node=node->first_node(); node; node=node->next_sibling())
			nodelist.push_back(node);
	}
}

unsigned long XMLfile::query(std::list<Node>& nodeselection, const std::string& querystring, Node startnode) const
{
	if (querystring.empty())
		// nothing selected
		return 0;
	// initialize startnode with current node of xmlfile object
	//const t_XMLnode* node(m_currentnode.m_xmlnode);
	// only ELEMENT_Node nodes will be handled => t_XMLelement* is the better choice (and doesn't need casts afterwards)...
	const t_XMLelement* node(static_cast<const t_XMLelement*>(m_currentnode.m_xmlnode));
	std::string nodepath(m_currentnode.nodepath());
	if(startnode.m_xmlnode)
	{ // if startnode is given use this one
		if(startnode.type()==Node::ELEMENT_Node)
		{ // but only if it's an ELEMENT_Node node
			node=static_cast<const t_XMLelement*>(startnode.m_xmlnode);
			nodepath=startnode.nodepath();
		}
		else
			std::cerr << "XMLfile::query: invalid startnode type " << startnode.type() << std::endl;
	}
	size_t pos=0;
	size_t tokenpos;
	if(querystring.find("/")==0)
	{ // querystr starts with "/" => absolute path: set startnode to document
		node=&m_xmldoc;
		nodepath=querystring[pos];
		++pos;
		if(querystring == "/") {
			nodeselection.push_back(Node(node,nodepath));
			return 1;
		}
	} else if(querystring.find("..")==0) {
		if(node->parent())
		{ // there exists a parent node
			node=node->parent();
			tokenpos=nodepath.find_last_of("/");
			if(tokenpos!=std::string::npos) nodepath=nodepath.substr(0,tokenpos+1);
			pos+=2;
			if(querystring[pos]=='/') ++pos;
		} else return 0;
	} else if(querystring.find(".")==0) {
		++pos;
		if(querystring[pos]=='/') ++pos;
	}
	unsigned long foundnodes=0;
	size_t pos2=querystring.find_first_of("/",pos+1);
	if(pos2==std::string::npos) pos2=querystring.size();
	std::string nodequery=querystring.substr(pos,pos2-pos);
	// nodequery is one chunk in the path separated by "/"
	//  assumed to be "elename" or "elename[condition]" for inner nodes and
	//  additionally "elename@attrname", "elename[condition]@attrname" for leaf nodes
	std::string elename, condition, attrname;
	tokenpos=nodequery.find_first_of("[@");
	unsigned long nodecondpos=0;
	std::string condattrname,condattrval;
	std::stringstream ss;
	if(tokenpos!=std::string::npos)
	{ // elename with condition and/or attribute
		elename=nodequery.substr(0,tokenpos);
		// for condition&attribute check condition first
		tokenpos=nodequery.find_first_of("[",tokenpos);
		if(tokenpos!=std::string::npos)
		{ // elename[condition]
			++tokenpos;
			size_t tokenpos2=nodequery.find_first_of("]",tokenpos);
			if(tokenpos2!=std::string::npos)
			{
				condition=nodequery.substr(tokenpos,tokenpos2-tokenpos);
				tokenpos=tokenpos2+1;
			} else {
				std::cerr << "ERROR: missing ] after " << elename << " within query string " << querystring << std::endl;
			}
			if(!condition.empty())
			{ // parse condition
				if(condition[0]=='@')
				{	// condition for attribute present
					tokenpos2=condition.find_first_of("=");
					condattrname=condition.substr(1,tokenpos2-1);
					condattrval=condition.substr(tokenpos2+1);
					tokenpos2=condattrval.find_first_of("\'\"");
					if(tokenpos2!=std::string::npos) condattrval=condattrval.substr(tokenpos2+1);
					tokenpos2=condattrval.find_last_of("\'\"");
					if(tokenpos2!=std::string::npos) condattrval=condattrval.substr(0,tokenpos2);
				} else { // node position specification is assumed as condition
					ss << condition;
					ss >> nodecondpos;
					ss.str(""); ss.clear();
				}
			}
		}
		else
		{
			tokenpos=0;
		}
		tokenpos=nodequery.find_first_of("@",tokenpos);
		if(tokenpos!=std::string::npos)
		{ // @attrname
			attrname=nodequery.substr(tokenpos+1);
		} else { // no attribute specification found
			attrname.clear();
		}
	} else { // just elename
		elename=nodequery;
		tokenpos=0;
	}
	// now nodequery should be separated into elename, condition, attrname

	std::string nodepath0(nodepath);
	// remove trailing /
	if(nodepath[nodepath.size()-1]=='/') nodepath0=nodepath.substr(0,nodepath.size()-1);
	t_XMLattribute* attr=NULL;
	unsigned long elepos=0;
	const t_XMLelement* ele;
	if(elename=="*")
	{ // wildcard *
		std::map<std::string,unsigned long> elenames;
		std::map<std::string,unsigned long>::iterator poselenames;
		for(ele=node->first_node(); ele; ele=ele->next_sibling())
		{
			poselenames=elenames.find(ele->name());
			if(poselenames==elenames.end())
				elepos=1;
			else
				elepos=poselenames->second+1;
			elenames[ele->name()]=elepos;
			ss << nodepath0 << "/" << ele->name() << "[" << elepos << "]";
			nodepath=ss.str();
			ss.str(""); ss.clear();
			if(nodecondpos)
			{	// check for node position condition
				if(elepos>nodecondpos) break;
				if(elepos!=nodecondpos) continue;
			}
			if(!condattrname.empty())
			{	// check for node attribute value
				attr=ele->first_attribute(condattrname.c_str());
				if(!attr || attr->value()!=condattrval) continue;
			}
			if(pos2<querystring.size()-1)
				foundnodes+=query(nodeselection,querystring.substr(pos2+1).c_str(),Node(ele,nodepath));
			else {
				nodeselection.push_back(Node(ele,nodepath));
				++foundnodes;
			}
		}
	} else { // elename given (or empty)
		const t_XMLelement* firstnode=node;
		if(!elename.empty())
			firstnode=firstnode->first_node(elename.c_str());
		for(ele=firstnode; ele; elename.empty() ? ele=NULL : ele=ele->next_sibling(elename.c_str()))
		{ // loop over all elements with given name
			++elepos;
			if(elename.empty())
			{
				nodepath=nodepath0;
			} else {
				ss << nodepath0 << "/" << ele->name() << "[" << elepos << "]";
				nodepath=ss.str();
				ss.str(""); ss.clear();
			}
			if(nodecondpos)
			{	// check for node position condition
				if(elepos>nodecondpos) break;
				if(elepos!=nodecondpos) continue;
			}
			if(!condattrname.empty())
			{	// check for node attribute value
				attr=ele->first_attribute(condattrname.c_str());
				if(!attr || attr->value()!=condattrval) continue;
			}
			if(pos2<querystring.size()-1)
			{ // node is inner node => recursion
				foundnodes+=query(nodeselection,querystring.substr(pos2+1).c_str(),Node(ele,nodepath));
			} else { // node found
				if (! attrname.empty())
				{ // search for attribute node
					attr=ele->first_attribute(attrname.c_str());
					nodepath.append("@");
					nodepath.append(attrname);
					nodeselection.push_back(Node(attr,nodepath));
				} else { // found element node
					nodeselection.push_back(Node(ele,nodepath));
				}
				++foundnodes;
			}
		}
	}
	return foundnodes;
}



void XMLfile::insertcloneelement(const t_XMLelement* src, t_XMLelement* dest_after)
{
	std::stringstream ss;
	ss << *src;
	char* src_content = m_xmldoc.allocate_string(ss.str().c_str());
	xml_document<> xmldoc;
	xmldoc.parse<0>(src_content);
	xml_node<>* cloneele=m_xmldoc.clone_node(xmldoc.first_node());
	xml_node<>* dest_parent=dest_after->parent();
	dest_parent->insert_node(dest_after,cloneele);
}


#ifdef ENABLE_MPI
void XMLfile::setMPIdefaults(int mpirootrank, MPI_Comm mpicomm)
{
	m_mpi_rootrank=mpirootrank;
	m_mpi_comm=mpicomm;

	int mpimyrank;
	MPI_CHECK( MPI_Comm_rank(m_mpi_comm, &mpimyrank) );
	m_mpi_myrank=mpimyrank;
	//m_mpi_isroot=mpirootrank==mpimyrank;
}


bool XMLfile::distributeXMLstring()
{
	std::string xmlstring;
	if(m_mpi_myrank==m_mpi_rootrank) xmlstring=std::string(*this);

	int ierror;
	int len=xmlstring.size();
	MPI_CHECK( ierror = MPI_Bcast(&len, 1, MPI_INT, m_mpi_rootrank, m_mpi_comm) );
	bool status = (ierror==MPI_SUCCESS) && (len>0);
	char* buffer=const_cast<char*>(xmlstring.c_str());
	if(m_mpi_myrank!=m_mpi_rootrank)
	{
		clear();
		buffer=m_xmldoc.allocate_string(NULL,len+1);
		buffer[len]=0;
	}
	MPI_CHECK( MPI_Bcast(buffer, len, MPI_CHAR, m_mpi_rootrank, m_mpi_comm) );
	if(m_mpi_myrank!=m_mpi_rootrank)
	{
		m_xmldoc.parse<0>(buffer);
		m_currentnode=Node(&m_xmldoc);
	}
	else
		status=false;

	return status;
}
#endif


// XMLfile::Node================================================================
// public methods

bool XMLfile::Node::isLeafNode() const
{
	if(!m_xmlnode) return false;
	if(type()!=ELEMENT_Node) return false;
	if(static_cast<const rapidxml::xml_node<>*>(m_xmlnode)->first_node()) return false;
	return true;
}

template<typename T> bool XMLfile::Node::getValue(T& value) const
{
	if(m_xmlnode)
	{
		std::stringstream ss(m_xmlnode->value());
		// Check if input has correct sign
		if (std::is_unsigned_v<T>) {
			if (ss.str().find_first_of("-") != std::string::npos) {
				std::cerr << "ERROR parsing \"" << ss.str() << "\" to data type " << typeid(T).name() << " from tag \"<" << name() << ">\" in xml file" << std::endl;
				std::cerr << "The tag contains a negative value but an unsigned value was expected." << std::endl;
				Simulation::exit(1);
			}
		}
		ss >> value;
		// Check if the entire string was consumed
		if (!ss.eof() || ss.fail()) {
			std::cerr << "ERROR parsing all chars of \"" << ss.str() << "\" from tag \"<" << name() << ">\" in xml file" << std::endl;
			std::cerr << "This might be the result of using a float while an integer is expected." << std::endl;
			Simulation::exit(1);
		}
		return true;
	}
	else
		return false;
}

template<> bool XMLfile::Node::getValue<std::string>(std::string& value) const
{
	if(m_xmlnode)
	{
		value=std::string(m_xmlnode->value());
		return true;
	}
	else
		return false;
}

template<> bool XMLfile::Node::getValue<bool>(bool& value) const
{
	std::string v;
	bool found=getValue(v);
	if (found) {
		// Remove white spaces
		v.erase(std::remove_if(v.begin(), v.end(), ::isspace), v.end());
		// Convert to upper case letters
		std::transform(v.begin(), v.end(), v.begin(), ::toupper);

		if (v == "TRUE" || v == "YES" || v == "ON") {
			value = true;
		} else if (v == "FALSE" || v == "NO" || v == "OFF") {
			value = false;
		} else {
			std::cerr << "ERROR parsing \"" << v << "\" to boolean from tag \"<" << name() << ">\" in xml file."
				<< "Valid values are: true, false, yes, no, on, off. " << std::endl;
			Simulation::exit(1);
		}
	}
	return found;
}

std::string XMLfile::Node::value_string(std::string defaultvalue) const
{
	std::string value(defaultvalue);
	if(m_xmlnode)
		getValue(value);
	else
		std::cerr << "XMLfile::Node::value_string: invalid node" << std::endl;
	return value;
}

int XMLfile::Node::value_int(int defaultvalue) const
{
	int value=defaultvalue;
	if(m_xmlnode)
		getValue(value);
	else
		std::cerr << "XMLfile::Node::value_int: invalid node" << std::endl;
	return value;
}

long XMLfile::Node::value_long(long defaultvalue) const
{
	long value=defaultvalue;
	if(m_xmlnode)
		getValue(value);
	else
		std::cerr << "XMLfile::Node::value_long: invalid node" << std::endl;
	return value;
}

float XMLfile::Node::value_float(float defaultvalue) const
{
	float value=defaultvalue;
	if(m_xmlnode)
		getValue(value);
	else
		std::cerr << "XMLfile::Node::value_float: invalid node" << std::endl;
	return value;
}

double XMLfile::Node::value_double(double defaultvalue) const
{
	double value=defaultvalue;
	if(m_xmlnode)
		getValue(value);
	else
		std::cerr << "XMLfile::Node::value_double: invalid node" << std::endl;
	return value;
}

bool XMLfile::Node::value_bool(bool defaultvalue) const
{
	bool value=defaultvalue;
	if(m_xmlnode)
		getValue(value);
	else
		std::cerr << "XMLfile::Node::value_bool: invalid node" << std::endl;
	return value;
}

void XMLfile::Node::printXML(std::ostream& ostrm) const
{
	if(*this)
	{
		if(m_type==ATTRIBUTE_Node)
		{
			std::cerr << "XMLfile::Node::printstrm: no XML representation of an attribute node present" << std::endl;
		} else {
			ostrm << *(static_cast<const t_XMLelement*>(m_xmlnode));
		}
	} else {
		std::cerr << "XMLfile::Node::printstrm: invalid node" << std::endl;
	}
}

void XMLfile::Node::print(std::ostream& ostrm) const
{
	if(*this)
	{
		ostrm << m_nodepath << ":\t" << m_xmlnode->name() << "=" << m_xmlnode->value() << std::endl;
	} else {
		ostrm << "invalid node" << std::endl;
	}
}

XMLfile::Node::Node(const t_XMLelement* xmlelement, std::string nodepath)
	: m_xmlnode(xmlelement), m_nodepath(nodepath), m_type(Unknown_Node)
{
	if(!xmlelement)
		m_type=Invalid_Node;
	else
		switch (xmlelement->type())
		{
			case node_document:
				m_type=DOCUMENT_Node;
				break;
			case node_element:
				m_type=ELEMENT_Node;
				break;
			case node_data:
				m_type=TEXT_Node;
				break;
			case node_cdata:
				m_type=CDATA_SECTION_Node;
				break;
			case node_comment:
				m_type=COMMENT_Node;
				break;
			case node_declaration:
				m_type=NOTATION_Node;
				break;
			case node_doctype:
				m_type=DOCUMENT_TYPE_Node;
				break;
			case node_pi:
				m_type=PROCESSING_INSTRUCTION_Node;
				break;
		}
}


// XMLfile::Query===============================================================
// public methods

XMLfile::Query& XMLfile::Query::operator =(const XMLfile::Query& q)
{
	xmlfile_unregister();
	m_xmlfile=q.m_xmlfile;
	m_nodes=q.m_nodes;
	xmlfile_register();
	return *this;
}

template<typename T> unsigned long XMLfile::Query::getNodeValue(T& value) const
{
	if(!empty()) front().getValue(value);
	return card();
}
template unsigned long XMLfile::Query::getNodeValue(int& value) const;
template unsigned long XMLfile::Query::getNodeValue(long& value) const;
template unsigned long XMLfile::Query::getNodeValue(unsigned int& value) const;
template unsigned long XMLfile::Query::getNodeValue(unsigned long& value) const;
template unsigned long XMLfile::Query::getNodeValue(unsigned long long& value) const;
template unsigned long XMLfile::Query::getNodeValue(float& value) const;
template unsigned long XMLfile::Query::getNodeValue(double& value) const;
template unsigned long XMLfile::Query::getNodeValue(bool& value) const;
template unsigned long XMLfile::Query::getNodeValue(std::string& value) const;

void XMLfile::Query::printXML(std::ostream& ostrm) const
{
	for(std::vector<Node>::const_iterator pos=m_nodes.begin();pos!=m_nodes.end();++pos)
		(*pos).printXML(ostrm);
}

void XMLfile::Query::print(std::ostream& ostrm) const
{
	for(std::vector<Node>::const_iterator pos=m_nodes.begin();pos!=m_nodes.end();++pos)
		(*pos).print(ostrm);
}
