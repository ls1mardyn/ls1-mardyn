#ifndef XMLREADER_MAIN_H_
#define XMLREADER_MAIN_H_

#include "External/tinyxpath/xpath_static.h"
#include "md_io/InputBase.h"
#include <string>
#include <iostream>

namespace md_io{
   class XMLReader_main;
}
using namespace std;

class md_io::XMLReader_main {
 public:
   XMLReader_main();
   ~XMLReader_main();
   TiXmlDocument XMLReader_main_get_doc(string filename);
   TiXmlDocument XMLReader_main_get_doc_chara(const char *xmldoc);

   string Eval_str(TiXmlDocument * XMLdoc_p, string expr);
   double Eval_d(TiXmlDocument * XMLdoc_p, string expr);
   int Eval_i(TiXmlDocument * XMLdoc_p, string expr);
   unsigned int Eval_ui(TiXmlDocument * XMLdoc_p, string expr);
   unsigned long Eval_ul(TiXmlDocument * XMLdoc_p, string expr);

   bool is_sourced(TiXmlDocument * XMLdoc_p, string expr);

   string merge(TiXmlDocument * XMLdoc_p, std::string pathPrefix);

 private:
   string _filename;
   TiXmlNode * traverse_descendents(TiXmlElement * element, std::string pathPrefix);

   TiXmlDocument * XMLdoc_g;

};

#endif /*XMLREADER_MAIN_H_*/
