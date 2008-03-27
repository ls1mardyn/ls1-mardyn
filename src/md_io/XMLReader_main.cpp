#include <string>
#include "md_io/XMLReader_main.h"

using namespace std;

md_io::XMLReader_main::XMLReader_main() {}

md_io::XMLReader_main::~XMLReader_main() {}

TiXmlDocument md_io::XMLReader_main::XMLReader_main_get_doc(string filename) {
   
   TiXmlDocument * XDp_doc;

   XDp_doc = new TiXmlDocument;
   if (! XDp_doc -> LoadFile (filename.c_str()))
   {    
      cout << "Failed to load \"" << filename << "\". No such file." << endl;
      std::cout << "error description: " << XDp_doc -> ErrorDesc () << std::endl;
      exit(1);
   }

   return *XDp_doc;
}

string md_io::XMLReader_main::Eval_str(TiXmlDocument *XMLdoc_p, string expr)
{
   TiXmlElement * XNp_root = XMLdoc_p -> RootElement ();
   TIXML_STRING res  = TinyXPath::S_xpath_string (XNp_root, expr.c_str());
   return res.c_str();
}

double md_io::XMLReader_main::Eval_d(TiXmlDocument *XMLdoc_p, string expr)
{
   TiXmlElement * XNp_root = XMLdoc_p -> RootElement ();
   TIXML_STRING res  = TinyXPath::S_xpath_string (XNp_root, expr.c_str());
   return strtod(res.c_str(),NULL);
}

int md_io::XMLReader_main::Eval_i(TiXmlDocument *XMLdoc_p, string expr)
{
   TiXmlElement * XNp_root = XMLdoc_p -> RootElement ();
   TIXML_STRING res  = TinyXPath::S_xpath_string (XNp_root, expr.c_str());
   return atoi(res.c_str());
}

unsigned int md_io::XMLReader_main::Eval_ui(TiXmlDocument *XMLdoc_p, string expr)
{
   TiXmlElement * XNp_root = XMLdoc_p -> RootElement ();
   TIXML_STRING res  = TinyXPath::S_xpath_string (XNp_root, expr.c_str());
   return atoi(res.c_str());
}

unsigned long md_io::XMLReader_main::Eval_ul(TiXmlDocument *XMLdoc_p, string expr)
{
   TiXmlElement * XNp_root = XMLdoc_p -> RootElement ();
   TIXML_STRING res  = TinyXPath::S_xpath_string (XNp_root, expr.c_str());
   return strtoul(res.c_str(),NULL,0);
}

bool md_io::XMLReader_main::is_sourced(TiXmlDocument * XMLdoc_p, string expr)
{
   TiXmlElement * XNp_root = XMLdoc_p -> RootElement ();
   string test_expr = (string)"count("+expr+"/@source)";
   TIXML_STRING res  = TinyXPath::S_xpath_string (XNp_root, test_expr.c_str());
   if (res == "1") {
	   return true;
   } else {
	   return false;
   }
}

TiXmlNode * md_io::XMLReader_main::traverse_descendents(TiXmlElement * element,
    std::string pathPrefix)
{
  // Recursively check all the descendents of the root node if it has
  // a source attribute. Each such node is replaced by the XML tree it
  // references to.
  TiXmlNode *node = element->IterateChildren(NULL);
  while (node)
  {
    if (node->Type() == 1)
    {
      const char * file = (node->ToElement())->Attribute("source");
      if (file != NULL)
      {
        string sfile(file);
        sfile = pathPrefix + sfile;
        // verify xml extension
        string::size_type loc1, loc2, loc3;
        loc1 = sfile.find(".xml", 0);
        loc2 = sfile.find(".XML", 0);
        loc3 = sfile.find(".Xml", 0);
        if ( (loc1 != string::npos ) || (loc2 != string::npos ) || (loc3
            != string::npos ))
        {
          // load the XML file which should replace the node
          TiXmlDocument new_source = XMLReader_main_get_doc(sfile);
          TiXmlDocument *new_source_p = &new_source;
          TiXmlElement *new_source_root = new_source_p->RootElement();
          // insert a proper attribute to the root node so that it
          // can be identified as a sourced node
          new_source_root->SetAttribute("format", "sourced");
          TiXmlNode *node_copy = new_source_root->Clone();

          // replace the node
          TiXmlNode *parent_node = node->Parent();
          TiXmlNode * return_node = parent_node->ReplaceChild(node, *node_copy);
          // check whether replacement was successful
          if (return_node == NULL)
          {
            cout << "Error while sourcing XML node!" << endl;
            exit(1);
          }
        }
      }
    }
    if ((node->NoChildren()==false) && (node->FirstChild()->Type() == 1))
    {
      traverse_descendents(node->ToElement(), pathPrefix);
    }
    node = element->IterateChildren(node);
  }
}

void md_io::XMLReader_main::merge(TiXmlDocument * XMLdoc_p, std::string pathPrefix)
{
   XMLdoc_g = XMLdoc_p;
   TiXmlElement *XNp_root = XMLdoc_g->RootElement();
  
   TiXmlNode *node = traverse_descendents(XNp_root, pathPrefix);
   
   XMLdoc_g->SaveFile("_temp.xml");
}
