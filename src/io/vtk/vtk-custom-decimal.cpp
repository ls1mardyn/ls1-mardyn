#include "xml-schema.h"

#include <limits>
#include <locale>
#include <sstream>

#include <xsd/cxx/ro-string.hxx>
#include <xsd/cxx/zc-istream.hxx>

using namespace std;

// Parsing
//
namespace xsd
{
  namespace cxx
  {
    namespace tree
    {
      long double
      traits<long double, char, schema_type::decimal>::create (const std::string& s,
              const xercesc::DOMElement*,
              flags,
              type*)
      {
      // This type cannot have whitespaces in its values. As result we
      // don't need to waste time collapsing whitespaces. All we need to
      // do is trim the string representation which can be done without
      // copying.
      //
      ro_string<char> tmp (s);
      trim (tmp);

      zc_istream<char> is (tmp);
      is.imbue (locale::classic ());

      long double t;
      is >> t;

      return t;
      }
    }
  }
}

//Serialization
//

namespace XERCES_CPP_NAMESPACE
{
  void
  operator<< (xercesc::DOMElement& e, const xml_schema::as_decimal& d)
  {
    ostringstream os;
    os.imbue (locale::classic ());
    os.precision (std::numeric_limits<long double>::digits10);

    os << std::fixed << d.x;

    // Remove the trailing zeros and the decimal point if necessary.
    //
    std::string s (os.str ());
    std::string::size_type size (s.size ()), n (size);

    for (; n > 0 && s[n - 1] == '0'; --n);

    if (n > 0 && s[n - 1] == '.')
      --n;

    if (n != size)
      s.resize (n);

    e << s;
  }

  void
  operator<< (xercesc::DOMAttr& a, const xml_schema::as_decimal& d)
  {
    ostringstream os;
    os.imbue (locale::classic ());
    os.precision (std::numeric_limits<long double>::digits10);

    os << std::fixed << d.x;

    // Remove the trailing zeros and the decimal point if necessary.
    //
    std::string s (os.str ());
    std::string::size_type size (s.size ()), n (size);

    for (; n > 0 && s[n - 1] == '0'; --n);

    if (n > 0 && s[n - 1] == '.')
      --n;

    if (n != size)
      s.resize (n);

    a << s;

  }
}

namespace xsd
{
  namespace cxx
  {
    namespace tree
    {
      void
      operator<< (xml_schema::list_stream& ls, const xml_schema::as_decimal& d)
      {
        
        ostringstream os;
        os.imbue (locale::classic ());
        os.precision (std::numeric_limits<long double>::digits10);

        os << std::fixed << d.x;
      
        // Remove the trailing zeros and the decimal point if necessary.
        //
        std::string s (os.str ());
        std::string::size_type size (s.size ()), n (size);

        for (; n > 0 && s[n - 1] == '0'; --n);

        if (n > 0 && s[n - 1] == '.')
          --n;

        if (n != size)
          s.resize (n);
	
	      ls.os_ << s;
      }
    }
  }
}
