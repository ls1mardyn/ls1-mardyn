// This code is a stub!

#include <string>
#include <iostream>
#include <fstream>
#include <time.h>
#include <climits>
#include "../src/External/tinyxpath/xpath_static.h"


int main(int argc, char *argv[], char *env[])
{

  if (argc < 4)
  {
    std::cout
        << "Syntax: moldy2mardyn <inputfile.dat> <inputfile.inp> <outputfile.xml>"
        << std::endl;
    exit(1);
  }
  // define the output filename
  std::string outfile_name = "new-";
  outfile_name.append(argv[2]);

  TiXmlDocument doc;
  TiXmlDeclaration * decl = new TiXmlDeclaration( "1.0", "UTF-8", "" );
  TiXmlElement * root = new TiXmlElement( "mardyncfg" );
  TiXmlElement * header = new TiXmlElement( "header" );
  TiXmlElement * version = new TiXmlElement( "version" );

  // apply the current date instead because we supposedly up to date
  time_t rawtime;
  struct tm * timeinfo;
  char datestr [10];
  time( &rawtime);
  timeinfo = localtime( &rawtime);
  strftime(datestr, 10, "%Y%m%d", timeinfo);
  TiXmlText * version_text = new TiXmlText( datestr );

  TiXmlElement * required_plugins = new TiXmlElement( "required-plugins" );
  TiXmlElement * experiment = new TiXmlElement( "experiment" );
  TiXmlElement * timestep_length = new TiXmlElement( "timestep-length" );
  TiXmlElement * cutoff_radius = new TiXmlElement( "cutoff-radius" );
  TiXmlElement * temperature = new TiXmlElement( "temperature" );
  TiXmlElement * current_time = new TiXmlElement( "current-time" );
  TiXmlElement * length = new TiXmlElement( "length" );
  TiXmlElement * lengthx = new TiXmlElement( "x" );
  TiXmlElement * lengthy = new TiXmlElement( "y" );
  TiXmlElement * lengthz = new TiXmlElement( "z" );
  TiXmlElement * phase_space = new TiXmlElement( "phase-space" );
  TiXmlElement * components = new TiXmlElement( "components" );
  TiXmlElement * components_data = new TiXmlElement( "data" );
  TiXmlElement * data_structure = new TiXmlElement( "data-structure" );

  root->LinkEndChild(header);
  header->LinkEndChild(version);
  version->LinkEndChild(version_text);
  header->LinkEndChild(required_plugins);
  root->LinkEndChild(experiment);
  experiment->LinkEndChild(timestep_length);
  experiment->LinkEndChild(cutoff_radius);
  experiment->LinkEndChild(temperature);
  experiment->LinkEndChild(current_time);
  experiment->LinkEndChild(length);
  length->LinkEndChild(lengthx);
  length->LinkEndChild(lengthy);
  length->LinkEndChild(lengthz);
  experiment->LinkEndChild(phase_space);
  experiment->LinkEndChild(components);
  components->LinkEndChild(components_data);
  experiment->LinkEndChild(data_structure);

  phase_space->SetAttribute("source", outfile_name.c_str() );
  phase_space->SetAttribute("format", "ASCII");

  std::string token;
  fstream inputfstream;

  // Parse the .inp file
  inputfstream.open(argv[2]);
  if (!inputfstream.is_open())
  {
    std::cout << "Error opening file " << argv[2] << " for reading.";
    exit(1);
  }

  while (inputfstream)
  {
    token.clear();
    inputfstream >> token;
    if (token.substr(0, 1)=="#")
    {
      inputfstream.ignore(INT_MAX, '\n');
    } else if (token.substr(0, 1)=="")
    {
      inputfstream.ignore(INT_MAX, '\n');
    } else if (token == "MOLDY")
    {
      inputfstream >> token;
    } else if (token == "t")
    {
      inputfstream >> token;
      TiXmlText * current_time_text = new TiXmlText( token.c_str() );
      current_time->LinkEndChild(current_time_text);
    } else if (token == "dt")
    {
      inputfstream >> token;
      TiXmlText * timestep_length_text = new TiXmlText( token.c_str() );
      timestep_length->LinkEndChild(timestep_length_text);
    } else if (token == "T")
    {
      inputfstream >> token;
      TiXmlText * temperature_text = new TiXmlText( token.c_str() );
      temperature->LinkEndChild(temperature_text);
      temperature->SetAttribute("unit", "kelvin");
    } else if (token == "L")
    {
      inputfstream >> token;
      TiXmlText * lengthx_text = new TiXmlText( token.c_str() );
      inputfstream >> token;
      TiXmlText * lengthy_text = new TiXmlText( token.c_str() );
      inputfstream >> token;
      TiXmlText * lengthz_text = new TiXmlText( token.c_str() );
      lengthx->LinkEndChild(lengthx_text);
      lengthy->LinkEndChild(lengthy_text);
      lengthz->LinkEndChild(lengthz_text);
    } else if (token == "rc")
    {
      inputfstream >> token;
      TiXmlText * cutoff_radius_text = new TiXmlText( token.c_str() );
      cutoff_radius->LinkEndChild(cutoff_radius_text);
    } else if (token == "C")
    {
      inputfstream >> token;
      int numcomponents = atoi(token.c_str());

      TiXmlElement * dcomponents = new TiXmlElement( "components" );
      dcomponents->SetAttribute("amount", numcomponents);
      components_data->LinkEndChild(dcomponents);

      std::string x, y, z, m, sigma, eps, xi, eta;
      unsigned int i, j;

      for (i=0; i<numcomponents; ++i)
      {
        unsigned int numljcenters=0;
        unsigned int numdipoles=0;
        unsigned int numquadrupoles=0;
        inputfstream >> numljcenters >> numdipoles >> numquadrupoles;
        TiXmlElement * comp = new TiXmlElement( "comp" );
        comp->SetAttribute("id", i+1);
        dcomponents->LinkEndChild(comp);

        for (j=0; j<numljcenters; ++j)
        {
          TiXmlElement * ljcenter = new TiXmlElement( "ljcenter" );
          ljcenter->SetAttribute("id", j+1);
          comp->LinkEndChild(ljcenter);

          inputfstream >> x >> y >> z >> m >> eps >> sigma;

          TiXmlElement * ljx = new TiXmlElement( "x" );
          TiXmlElement * ljy = new TiXmlElement( "y" );
          TiXmlElement * ljz = new TiXmlElement( "z" );
          TiXmlElement * ljm = new TiXmlElement( "m" );
          TiXmlElement * ljeps = new TiXmlElement( "eps" );
          TiXmlElement * ljsigma = new TiXmlElement( "sigma" );

          TiXmlText * ljx_text = new TiXmlText( x.c_str() );
          TiXmlText * ljy_text = new TiXmlText( y.c_str() );
          TiXmlText * ljz_text = new TiXmlText( z.c_str() );
          TiXmlText * ljm_text = new TiXmlText( m.c_str() );
          TiXmlText * ljeps_text = new TiXmlText( eps.c_str() );
          TiXmlText * ljsigma_text = new TiXmlText( sigma.c_str() );

          ljx->LinkEndChild(ljx_text);
          ljy->LinkEndChild(ljy_text);
          ljz->LinkEndChild(ljz_text);
          ljm->LinkEndChild(ljm_text);
          ljeps->LinkEndChild(ljeps_text);
          ljsigma->LinkEndChild(ljsigma_text);

          ljcenter->LinkEndChild(ljx);
          ljcenter->LinkEndChild(ljy);
          ljcenter->LinkEndChild(ljz);
          ljcenter->LinkEndChild(ljm);
          ljcenter->LinkEndChild(ljeps);
          ljcenter->LinkEndChild(ljsigma);
        }
        for (j=0; j<numdipoles; ++j)
        {
          std::string eMyx, eMyy, eMyz, absMy;
          TiXmlElement * dipole = new TiXmlElement( "dipole" );
          dipole->SetAttribute("id", j+1);
          comp->LinkEndChild(dipole);

          inputfstream >> x >> y >> z >> eMyx >> eMyy >> eMyz >> absMy;

          TiXmlElement * dix = new TiXmlElement( "x" );
          TiXmlElement * diy = new TiXmlElement( "y" );
          TiXmlElement * diz = new TiXmlElement( "z" );
          TiXmlElement * dieMyx = new TiXmlElement( "eMyx" );
          TiXmlElement * dieMyy = new TiXmlElement( "eMyy" );
          TiXmlElement * dieMyz = new TiXmlElement( "eMyz" );
          TiXmlElement * diabsMy = new TiXmlElement( "absMy" );

          TiXmlText * dix_text = new TiXmlText( x.c_str() );
          TiXmlText * diy_text = new TiXmlText( y.c_str() );
          TiXmlText * diz_text = new TiXmlText( z.c_str() );
          TiXmlText * dieMyx_text = new TiXmlText( eMyx.c_str() );
          TiXmlText * dieMyy_text = new TiXmlText( eMyy.c_str() );
          TiXmlText * dieMyz_text = new TiXmlText( eMyz.c_str() );
          TiXmlText * diabsMy_text = new TiXmlText( absMy.c_str() );

          dix->LinkEndChild(dix_text);
          diy->LinkEndChild(diy_text);
          diz->LinkEndChild(diz_text);
          dieMyx->LinkEndChild(dieMyx_text);
          dieMyy->LinkEndChild(dieMyy_text);
          dieMyz->LinkEndChild(dieMyz_text);
          diabsMy->LinkEndChild(diabsMy_text);

          dipole->LinkEndChild(dix);
          dipole->LinkEndChild(diy);
          dipole->LinkEndChild(diz);
          dipole->LinkEndChild(dieMyx);
          dipole->LinkEndChild(dieMyy);
          dipole->LinkEndChild(dieMyz);
          dipole->LinkEndChild(diabsMy);
        }
        for (j=0; j<numquadrupoles; ++j)
        {
          std::string eQx, eQy, eQz, absQ;
          TiXmlElement * quadrupole = new TiXmlElement( "quadrupole" );
          quadrupole->SetAttribute("id", j+1);
          comp->LinkEndChild(quadrupole);

          inputfstream >> x >> y >> z >> eQx >> eQy >> eQz >> absQ;

          TiXmlElement * ljx = new TiXmlElement( "x" );
          TiXmlElement * ljy = new TiXmlElement( "y" );
          TiXmlElement * ljz = new TiXmlElement( "z" );
          TiXmlElement * ljeQx = new TiXmlElement( "eQx" );
          TiXmlElement * ljeQy = new TiXmlElement( "eQy" );
          TiXmlElement * ljeQz = new TiXmlElement( "eQz" );
          TiXmlElement * ljabsQ = new TiXmlElement( "absQ" );

          TiXmlText * ljx_text = new TiXmlText( x.c_str() );
          TiXmlText * ljy_text = new TiXmlText( y.c_str() );
          TiXmlText * ljz_text = new TiXmlText( z.c_str() );
          TiXmlText * ljeQx_text = new TiXmlText( eQx.c_str() );
          TiXmlText * ljeQy_text = new TiXmlText( eQy.c_str() );
          TiXmlText * ljeQz_text = new TiXmlText( eQz.c_str() );
          TiXmlText * ljabsQ_text = new TiXmlText( absQ.c_str() );

          ljx->LinkEndChild(ljx_text);
          ljy->LinkEndChild(ljy_text);
          ljz->LinkEndChild(ljz_text);
          ljeQx->LinkEndChild(ljeQx_text);
          ljeQy->LinkEndChild(ljeQy_text);
          ljeQz->LinkEndChild(ljeQz_text);
          ljabsQ->LinkEndChild(ljabsQ_text);

          quadrupole->LinkEndChild(ljx);
          quadrupole->LinkEndChild(ljy);
          quadrupole->LinkEndChild(ljz);
          quadrupole->LinkEndChild(ljeQx);
          quadrupole->LinkEndChild(ljeQy);
          quadrupole->LinkEndChild(ljeQz);
          quadrupole->LinkEndChild(ljabsQ);
        }

        std::string IDummy1, IDummy2, IDummy3;
        inputfstream >> IDummy1 >> IDummy2 >> IDummy3;

        TiXmlElement * dummy1 = new TiXmlElement( "dummy1" );
        TiXmlElement * dummy2 = new TiXmlElement( "dummy2" );
        TiXmlElement * dummy3 = new TiXmlElement( "dummy3" );

        TiXmlText * dummy1_text = new TiXmlText( IDummy1.c_str() );
        TiXmlText * dummy2_text = new TiXmlText( IDummy2.c_str() );
        TiXmlText * dummy3_text = new TiXmlText( IDummy3.c_str() );

        dummy1->LinkEndChild(dummy1_text);
        dummy2->LinkEndChild(dummy2_text);
        dummy3->LinkEndChild(dummy3_text);

        comp->LinkEndChild(dummy1);
        comp->LinkEndChild(dummy2);
        comp->LinkEndChild(dummy3);

      }

      TiXmlElement * mixcoeff = new TiXmlElement( "mixcoeff" );
      dcomponents->LinkEndChild(mixcoeff);

      for (i=0; i<numcomponents-1; ++i)
      {
        for (j=i+1; j<numcomponents; ++j)
        {
          inputfstream >> xi >> eta;
          TiXmlElement * mix_xi = new TiXmlElement( "xi" );
          TiXmlElement * mix_eta = new TiXmlElement( "eta" );
          mix_xi->SetAttribute("id", i); // can be ignored
          mix_eta->SetAttribute("id", i); // can be ignored
          TiXmlText * mix_xi_text = new TiXmlText( xi.c_str() );
          TiXmlText * mix_eta_text = new TiXmlText( eta.c_str() );
          mix_xi->LinkEndChild(mix_xi_text);
          mix_eta->LinkEndChild(mix_eta_text);
          mixcoeff->LinkEndChild(mix_xi);
          mixcoeff->LinkEndChild(mix_eta);
        }
      }

      inputfstream >> token;
      TiXmlElement * epsilon_rf = new TiXmlElement( "epsilon-rf" );
      TiXmlText * epsilon_rf_text = new TiXmlText( token.c_str() );
      epsilon_rf->LinkEndChild(epsilon_rf_text);
      dcomponents->LinkEndChild(epsilon_rf);

      std::cout << "Warning: falling back to linked cells data-structure" << std::endl;
      TiXmlElement * data_structure_type = new TiXmlElement( "linked-cells" );
      data_structure->LinkEndChild(data_structure_type);
      TiXmlText * linked_cells_text = new TiXmlText( "1" );
      data_structure_type->LinkEndChild(linked_cells_text);

      std::ofstream outfile;
      outfile.open(outfile_name.c_str() );
      if (outfile.is_open())
      {
        getline(inputfstream, token);
        inputfstream >> token;
        inputfstream >> token;
        outfile << "NumberOfMolecules\t" << token << std::endl;
        inputfstream >> token;
        outfile << "MoleculeFormat\t\t" << token;
        while (inputfstream)
        {
          getline(inputfstream, token);
          outfile << token << std::endl;
        }
        outfile.close();
      } else
      {
        std::cout << "Unable to open file " << outfile_name << " for writing.";
        exit(1);
      }
    } else
    {
      std::cerr << "Warning: Unknown token '" << token << "' at position "
          << inputfstream.tellg() << std::endl;
    }
  }

  inputfstream.close();

  std::cout << " * the phase space point file has been saved in '" << outfile_name
      << "'" << std::endl;

  // Write the resulting XML tree to the file specified
  doc.LinkEndChild(decl);
  doc.LinkEndChild(root);
  doc.SaveFile(argv[3]);
}
