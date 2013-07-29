/** \file Expression.cpp
  * \brief Expression (tree) with variables
  * \author Martin Bernreuther <bernreuther@hlrs.de>
  * license: GPL (http://www.gnu.de/documents/gpl.de.html)
*/

#include "Expression.h"

using namespace std;

Expression::Variable* Expression::VariableSet::addVariable(const string& name)
{
	string vgrpname,varname(name);
	size_t colonpos=name.find(":");
	if(colonpos!=string::npos)
	{
		vgrpname=name.substr(0,colonpos);
		varname=name.substr(colonpos+1);
	}
	if(!existVariableGroup(vgrpname))
	{ // create variable group
		//_vargroups[vgrpname]=VariableGroup(vgrpname);
		_vargroups.insert(pair<string,VariableGroup>(vgrpname,VariableGroup(vgrpname)));
	}
	//_variables[name]=Variable(varname,&_vargroups[vgrpname],val);
	_variables.insert(pair<string,Variable>(name,Variable(varname,&_vargroups[vgrpname])));
	return &_variables[name];
}

bool Expression::VariableSet::removeVariable(const string& name)
{
	Variable* var=getVariable(name);
	if(var)
	{
		const VariableGroup* vargrp=var->getVariableGroup();
		_variables.erase(var->getfullName());
		if(vargrp)
		{
			//vargrp->removeVariable(name);	// done by Variable destructor
			if(vargrp->countVariables()) _vargroups.erase(vargrp->getName());
		}
		return true;
	}
	else
		return false;
}



void Expression::Node::traverse(list<const Node*>& nodelist, enum Etraversetype traversetype) const
{
	switch(traversetype)
	{
		case traversetypePREFIX:
			nodelist.push_back(this);
			if(_children[0]) _children[0]->traverse(nodelist,traversetype);
			if(_children[1]) _children[1]->traverse(nodelist,traversetype);
			break;
		case traversetypeINFIX:
			if(_children[0]) _children[0]->traverse(nodelist,traversetype);
			nodelist.push_back(this);
			if(_children[1]) _children[1]->traverse(nodelist,traversetype);
			break;
		case traversetypePOSTFIX:
			if(_children[0]) _children[0]->traverse(nodelist,traversetype);
			if(_children[1]) _children[1]->traverse(nodelist,traversetype);
			nodelist.push_back(this);
			break;
	}
}

void Expression::Node::writeSubExpr(ostream& ostrm, enum Etraversetype traversetype, char sep) const
{
	switch(traversetype)
	{
		case traversetypePREFIX:
			write(ostrm);
			if(sep) ostrm << sep;
			if(_children[0]) _children[0]->writeSubExpr(ostrm,traversetype,sep);
			if(_children[1]) _children[1]->writeSubExpr(ostrm,traversetype,sep);
			break;
		case traversetypeINFIX:
			if(_children[0])
			{
				if(_children[0]->_priority>_priority) ostrm<<"(";
				_children[0]->writeSubExpr(ostrm,traversetype,sep);
				if(_children[0]->_priority>_priority) ostrm<<")";
			}
			write(ostrm);
			if(sep) ostrm << sep;
			if(_children[1])
			{
				if(_children[1]->_priority>_priority) ostrm<<"(";
				_children[1]->writeSubExpr(ostrm,traversetype,sep);
				if(_children[1]->_priority>_priority) ostrm<<")";
			}
			break;
		case traversetypePOSTFIX:
			if(_children[0]) _children[0]->writeSubExpr(ostrm,traversetype,sep);
			if(_children[1]) _children[1]->writeSubExpr(ostrm,traversetype,sep);
			write(ostrm);
			if(sep) ostrm << sep;
			break;
	}
}


Expression::NodeFunction1::NodeFunction1(const string& name, Node* child, Node* parent)
	: Node(NULL,child,parent,0)
{
	if(name.compare("abs")==0||name.compare("ABS")==0)
		_functype=functypeABS;
	else if(name.compare("floor")==0||name.compare("FLOOR")==0)
		_functype=functypeFLOOR;
	else if(name.compare("ceil")==0||name.compare("CEIL")==0)
		_functype=functypeCEIL;
	else if(name.compare("sqrt")==0||name.compare("SQRT")==0)
		_functype=functypeSQRT;
	else if(name.compare("log")==0||name.compare("LOG")==0)
		_functype=functypeLOG;
	else if(name.compare("log2")==0||name.compare("LOG2")==0)
		_functype=functypeLOG2;
	else if(name.compare("log10")==0||name.compare("LOG10")==0)
		_functype=functypeLOG10;
	else if(name.compare("sin")==0||name.compare("SIN")==0)
		_functype=functypeSIN;
	else if(name.compare("cos")==0||name.compare("COS")==0)
		_functype=functypeCOS;
	else if(name.compare("tan")==0||name.compare("TAN")==0)
		_functype=functypeTAN;
	else
		_functype=functypeNONE;
}

Expression::Tvaltype Expression::NodeFunction1::valueType() const
{
	if(_children[1])
	{
		switch(_functype)
		{
			case functypeNONE: return valtypeNONE;
			case functypeABS: return _children[1]->valueType();
			case functypeFLOOR:
			case functypeCEIL:
				return valtypeINT;
			case functypeSQRT:
			case functypeLOG:
			case functypeLOG2:
			case functypeLOG10:
			case functypeSIN:
			case functypeCOS:
			case functypeTAN:
			//default:
				return valtypeFLOAT;
			//default: return valtypeNONE;
		}
	}
	return valtypeNONE;
}

Expression::Value Expression::NodeFunction1::evaluate() const
{
	if(_children[1])
	{
		switch(_functype)
		{
			case functypeNONE: return Value();
			case functypeABS:
				if(_children[1]->isInt())
					return Value(labs(Tint(_children[1]->evaluate())));
				else
					return Value(fabs(Tfloat(_children[1]->evaluate())));
			case functypeFLOOR: return Value(Tint(floor(Tfloat(_children[1]->evaluate()))));
			case functypeCEIL:  return Value(Tint(ceil(Tfloat(_children[1]->evaluate()))));
			case functypeSQRT:  return Value(sqrt(Tfloat(_children[1]->evaluate())));
			case functypeLOG:   return Value(log(Tfloat(_children[1]->evaluate())));
			case functypeLOG2:  return Value(log2(Tfloat(_children[1]->evaluate())));
			case functypeLOG10: return Value(log10(Tfloat(_children[1]->evaluate())));
			case functypeSIN:   return Value(sin(Tfloat(_children[1]->evaluate())));
			case functypeCOS:   return Value(cos(Tfloat(_children[1]->evaluate())));
			case functypeTAN:   return Value(tan(Tfloat(_children[1]->evaluate())));
			//case functypeNONE:
			default:
				return Value();
		}
	}
	else return 0.;
}

void Expression::NodeFunction1::write(ostream& ostrm) const {
	switch(_functype)
	{
		case functypeNONE:  ostrm << "undef"; break;
		case functypeABS:   ostrm << "ABS"; break;
		case functypeFLOOR: ostrm << "FLOOR"; break;
		case functypeCEIL:  ostrm << "CEIL"; break;
		case functypeSQRT:  ostrm << "SQRT"; break;
		case functypeLOG:   ostrm << "LOG"; break;
		case functypeLOG2:  ostrm << "LOG2"; break;
		case functypeLOG10: ostrm << "LOG10"; break;
		case functypeSIN:   ostrm << "SIN"; break;
		case functypeCOS:   ostrm << "COS"; break;
		case functypeTAN:   ostrm << "TAN"; break;
	}
}


Expression::NodeOperation2::NodeOperation2(char op, Node* child1, Node* child2, Node* parent)
	: Node(child1,child2,parent), _operator(op)
{
	switch(op)
	{
		case '*':
		case '/':
			_priority=2;
			break;
		case '-':
			_priority=3;
			break;
		case '+':
			_priority=4;
			break;
	}
}

Expression::Value Expression::NodeOperation2::evaluate() const
{
	if(_children[0] && _children[1])
	{
		switch(_operator)
		{
			case '+':
				if(_children[0]->isInt()&&_children[1]->isInt())
					return Value(_children[0]->evaluateInt()+_children[1]->evaluateInt());
				else if (_children[0]->isInt()&&_children[1]->isFloat())
					return Value(_children[0]->evaluateInt()+_children[1]->evaluateFloat());
				else if (_children[0]->isFloat()&&_children[1]->isInt())
					return Value(_children[0]->evaluateFloat()+_children[1]->evaluateInt());
				else if (_children[0]->isFloat()&&_children[1]->isFloat())
					return Value(_children[0]->evaluateFloat()+_children[1]->evaluateFloat());
				break;
			case '-':
				if(_children[0]->isInt()&&_children[1]->isInt())
					return Value(_children[0]->evaluateInt()-_children[1]->evaluateInt());
				else if (_children[0]->isInt()&&_children[1]->isFloat())
					return Value(_children[0]->evaluateInt()-_children[1]->evaluateFloat());
				else if (_children[0]->isFloat()&&_children[1]->isInt())
					return Value(_children[0]->evaluateFloat()-_children[1]->evaluateInt());
				else if (_children[0]->isFloat()&&_children[1]->isFloat())
					return Value(_children[0]->evaluateFloat()-_children[1]->evaluateFloat());
				break;
			case '*':
				if(_children[0]->isInt()&&_children[1]->isInt())
					return Value(_children[0]->evaluateInt()*_children[1]->evaluateInt());
				else if (_children[0]->isInt()&&_children[1]->isFloat())
					return Value(_children[0]->evaluateInt()*_children[1]->evaluateFloat());
				else if (_children[0]->isFloat()&&_children[1]->isInt())
					return Value(_children[0]->evaluateFloat()*_children[1]->evaluateInt());
				else if (_children[0]->isFloat()&&_children[1]->isFloat())
					return Value(_children[0]->evaluateFloat()*_children[1]->evaluateFloat());
				break;
			case '/':
				if(_children[0]->isInt()&&_children[1]->isInt())
					return Value(_children[0]->evaluateInt()/_children[1]->evaluateInt());
				else if (_children[0]->isInt()&&_children[1]->isFloat())
					return Value(_children[0]->evaluateInt()/_children[1]->evaluateFloat());
				else if (_children[0]->isFloat()&&_children[1]->isInt())
					return Value(_children[0]->evaluateFloat()/_children[1]->evaluateInt());
				else if (_children[0]->isFloat()&&_children[1]->isFloat())
					return Value(_children[0]->evaluateFloat()/_children[1]->evaluateFloat());
				break;
		}
	}
	return Value();
}



void Expression::initializeRPN(const string& exprstr, bool genlabel)
{
	clear();
	stack<Node*> nodestack;
	// split expr to generate Nodes
	size_t startpos=0,endpos;
	while(startpos<exprstr.size()) {
		endpos=exprstr.find_first_of(" \t",startpos);
		if(endpos==string::npos)
		{
			endpos=exprstr.size()-1;
		} else if(endpos==startpos) {
			startpos=endpos+1;
			continue;
		} else {
			--endpos;
		}
		string token=exprstr.substr(startpos,endpos-startpos+1);
		size_t colonpos=token.find(":");
		if(token.size()==1 && token.find_first_of("+-*/")==0)
		{	// operator ................................................
			char op=token[0];
			if(nodestack.size()>=2)
			{
				Node* node1=nodestack.top(); nodestack.pop();
				Node* node2=nodestack.top(); nodestack.pop();
				Node* node0=new NodeOperation2(op,node2,node1);
				nodestack.push(node0);
				node1->setParent(node0);
				node2->setParent(node0);
			}
		} else if(token.find_first_not_of("0123456789.-E")==string::npos){
			// constant ................................................
			istringstream iss(token);
			if(token.find_first_not_of("0123456789-")==string::npos)
			{ // Tint
				Tint valInt=0;
				iss>>valInt;
				nodestack.push(new NodeConstant(valInt));
			}
			else
			{ // Tfloat
				Tfloat valFloat=0.;
				iss>>valFloat;
				nodestack.push(new NodeConstant(valFloat));
			}
			iss.str(string());
		} else if(colonpos!=string::npos) {
			// variable ................................................
			nodestack.push(new NodeVariable(_variableset->addVariable(token)));
		} else {
			// function ................................................
			//string funcname=token;
			Node* node1=nodestack.top(); nodestack.pop();
			Node* node0=new NodeFunction1(token,node1);
			nodestack.push(node0);
			node1->setParent(node0);
		}
		startpos=endpos+2;
	}
	if(!nodestack.empty()) _rootnode=nodestack.top();
	if(genlabel) genLabel();
}

