/** \file Expression.h
  * \brief Expression (tree) with variables
  * \author Martin Bernreuther <bernreuther@hlrs.de>
*/

#ifndef EXPRESSION_H
#define EXPRESSION_H

#include <map>
#include <set>
#include <stack>
#include <list>
#include <string>
#include <sstream>
#include <iostream>

#include <cmath>
#include <cstdlib>	// abs

// Expression ------------------------------------------------------------------------------
class Expression
{
	public:
		typedef long Tint;
		typedef double Tfloat;

		enum Evaltype {valtypeNONE, valtypeINT, valtypeFLOAT};
		typedef enum Evaltype Tvaltype;
		enum Etraversetype { traversetypePREFIX, traversetypeINFIX, traversetypePOSTFIX };
		typedef enum Etraversetype Ttraversetype;

	// Value -----------------------------------------------------------------------------------
	class Value
	{
		public:
			union TUvalue
			{
				Tint valInt;
				Tfloat valFloat;
			};
			typedef union TUvalue Tvalue;

			/// Constructor
			Value() : _type(valtypeNONE) { _value.valInt=0; }
			/// Constructor
			/**
				parameter	type	Evaltype	type (INT/FLOAT)
				parameter	Tvalue	Tfloat	value
			**/
			Value(enum Evaltype type, Tvalue value)
			     : _type(type), _value(value) {}
			/// Constructor
			/**
				parameter	valInt		integer value
			**/
			//Value(Tint valInt) { _type=valtypeINT; _value.valInt=valInt; }
			Value(int valInt) { _type=valtypeINT; _value.valInt=Tint(valInt); }
			Value(long valInt) { _type=valtypeINT; _value.valInt=Tint(valInt); }
			/// Constructor
			/**
				parameter		Tfloat	floating point real value
			**/
			//Value(Tfloat valFloat) { _type=valtypeFLOAT; _value.valFloat=valFloat; }
			Value(float valFloat) { _type=valtypeFLOAT; _value.valFloat=Tfloat(valFloat); }
			Value(double valFloat) { _type=valtypeFLOAT; _value.valFloat=Tfloat(valFloat); }

			/// get value type
			/**
				parameter		Evaltype
				retval	Tvaltype	value type
			**/
			Tvaltype getType() const { return _type; }
			/// get value
			/**
				retval		Tvalue	value
			**/
			Tvalue getValue() const { return _value; }
			/// is the value an integer value?
			/**
				retval		bool	result
			**/
			bool isInt() const { return getType()==valtypeINT; }
			/// is the value an floating point value?
			/**
				retval		bool	result
			**/
			bool isFloat() const { return getType()==valtypeFLOAT; }
			/// get value as floating point number
			/**
				retval		Tfloat	value
			**/
			Tfloat getValueFloat() const
			{
				switch(_type)
				{
					case valtypeFLOAT: return _value.valFloat;
					case valtypeINT: return Tfloat(_value.valInt);
					default: return 0.;
				}
			}
			operator Tfloat() const { return getValueFloat(); }
			/// get value as integer number
			/**
				retval		Tint	value
			**/
			Tint getValueInt() const
			{
				switch(_type)
				{
					case valtypeINT: return _value.valInt;
					case valtypeFLOAT: return Tint(_value.valFloat);
					default: return 0.;
				}
			}
			operator Tint() const { return getValueInt(); }
			/// write value to stream
			/**
				parameter		std::ostrm&	stream
			**/
			void write(std::ostream& ostrm=std::cout) const
			{
				if(isInt())
					ostrm << getValueInt();
				else if(isFloat())
					ostrm << getValueFloat();
			}
			operator std::string() const
			{
				std::ostringstream oss;
				write(oss);
				return oss.str();
			}
			bool operator==(Value const& v) const;
			bool operator>(Value const& v) const;
			bool operator>=(Value const& v) const;
			bool operator<(Value const& v) const;
			bool operator<=(Value const& v) const;
			const Value operator+(Value const& v) const;
			const Value operator-(Value const& v) const;
			const Value operator*(Value const& v) const;
			const Value operator/(Value const& v) const;

		private:
			Tvaltype _type;
			Tvalue _value;
	};
	// ----------------------------------------------------------------------------------- Value

	class Variable;
	// VariableGroup ---------------------------------------------------------------------------
	/**
		A VariableSet might contain Variable- and VariableGroup-Objects,
		whereas a Variable might belong to a VariableGroup, which is also
		reflected by the Variable name: <VariableGroup>:<Variable>.
		The VariableGroup-map acts like a hash table
	**/
	class VariableGroup
	{
		public:
			VariableGroup() : _name(std::string()) {}
			VariableGroup(const std::string& name) : _name(name) {}
			const std::string getName() const { return _name; }
			void addVariable(const Variable* var) { _variables.insert(var); }
			bool removeVariable(const Variable* var)
			{	if(_variables.count(var)) { _variables.erase(var); return true; } else return false;	}
			unsigned int countVariables() const { return _variables.size(); }
			operator unsigned int() const { return countVariables(); }
		private:
			std::string _name;
			std::set<const Variable*> _variables;
	};
	// --------------------------------------------------------------------------- VariableGroup

	// Variable --------------------------------------------------------------------------------
	class Variable
	{
		public:
			Variable() : _name(std::string()), _value(Value()), _vargrp(NULL) {}
			Variable(const std::string& name, VariableGroup* grp)
			            : _name(name), _value(Value()), _vargrp(grp) { if(_vargrp) _vargrp->addVariable(this); }
			template <class T> Variable(const std::string& name, VariableGroup* grp, T val)
			            : _name(name), _value(val), _vargrp(grp) { if(_vargrp) _vargrp->addVariable(this); }
			~Variable() { if(_vargrp) _vargrp->removeVariable(this); }
			/// variable name
			/**
				retval		std::string	variable name
			**/
			const std::string getName() const { return _name; }
			/// full variable name
			/**
				retval		std::string	group:variable name
			**/
			const std::string getfullName() const
			{
				if(_vargrp)
					return _vargrp->getName()+":"+_name;
				else
					return _name;
			}
			Tvaltype getType() const { return _value.getType(); }
			Value getValue() const { return _value; }
			Tfloat getValueFloat() const { return _value.getValueFloat(); }
			operator Tfloat() const { return getValueFloat(); }
			Tint getValueInt() const { return _value.getValueInt(); }
			operator Tint() const { return getValueInt(); }
			bool isInt() const { return _value.isInt(); }
			bool isFloat() const { return _value.isFloat(); }
			const VariableGroup* getVariableGroup() const { return _vargrp; }
			/// set value
			/**
				parameter		T	value
			**/
			template <class T> void setValue(T val) { _value=Value(val); }
			void write(std::ostream& ostrm=std::cout, bool prtval=false) const
			{
				ostrm << getfullName();
				if(prtval)
				{
					ostrm << "{=";
					_value.write(ostrm);
					ostrm << "}";
				}
			}
			operator std::string() const
			{
				std::ostringstream oss;
				write(oss);
				return oss.str();
			}

		private:
			std::string _name;
			Value _value;
			VariableGroup* _vargrp;
	};
	// -------------------------------------------------------------------------------- Variable

	// VariableSet -----------------------------------------------------------------------------
	class VariableSet
	{
		public:
			VariableSet() {};

			const std::set<const VariableGroup*> getVariableGroupNames() const
			{
				std::set<const VariableGroup*> vargroups;
				for (std::map<std::string,VariableGroup>::const_iterator it=_vargroups.begin(); it!=_vargroups.end(); ++it)
					vargroups.insert(&it->second);
				return vargroups;
			}
			unsigned int VariableGroupsCount() const { return _vargroups.size(); }
			bool existVariableGroup(const std::string& name) const { return _vargroups.count(name)>0; }
			const VariableGroup* getVariableGroup(const std::string& name) const
			{	if(existVariableGroup(name)) return &(_vargroups.find(name)->second); else return NULL;	}
			unsigned int VariableGroupVariablesCount(const std::string& name) const
			{
				if(existVariableGroup(name))
					//return _vargroups[name].countVariables();
					return (_vargroups.find(name)->second).countVariables();
				else
					return 0;

			}
			unsigned int VariablesCount() const { return _variables.size(); }
			bool existVariable(const std::string& name) const { return _variables.count(name)>0; }
			std::set<Variable*> getVariables()
			{
				std::set<Variable*> variables;
				for(std::map<std::string,Variable>::iterator it=_variables.begin(); it!=_variables.end(); ++it)
					variables.insert(&it->second);
				return variables;
			}
			Variable* getVariable(const std::string& name)
			{
				if(existVariable(name))
					return &_variables[name];
				else
					return NULL;
			}
			Variable* addVariable(const std::string& name);
			template <class T> bool setVariable(const std::string& name,T val=0)
			{
				Variable* var=getVariable(name);
				if(var)
				{ // variable already exists
					var->setValue(val);
					return false;
				}
				else
				{ // create variable
					var=addVariable(name);
					var->setValue(val);
					return true;
				}
			}
			//template <class T> bool setVariable(const char* name,T val=0) { return setVariable(std::string(name),val); }
			template <class T> bool setVariable(const std::string& vgrpname, const std::string& varname, T val=0)
			{	return setVariable(std::string(vgrpname+":"+varname),val);	}
			//template <class T> bool setVariable(const char* vgrpname, const char* varname, T val=0) { return setVariable(std::string(vgrpname),std::string(varname),val); }
			bool removeVariable(const std::string& name);

		private:
			std::map<std::string,Variable> _variables;
			std::map<std::string,VariableGroup> _vargroups;
	};
	// ----------------------------------------------------------------------------- VariableSet

	// Node and derivatives --------------------------------------------------------------------

	// Node ------------------------------------------------------------------------------------
	class Node
	{
		public:
			Node(Node* child1=NULL, Node* child2=NULL, Node* parent=NULL, short priority=0)
			        : _parent(parent), _priority(priority)
			{
				_children[0]=child1;
				_children[1]=child2;
			}
			virtual ~Node()
			{
				if(_children[0]) delete(_children[0]);
				if(_children[1]) delete(_children[1]);
			}
			void setChild1(Node* child1) { _children[0]=child1; }
			Node* getChild1() const { return _children[0]; }
			void setChild2(Node* child2) { _children[1]=child2; }
			Node* getChild2() const { return _children[1]; }
			void setParent(Node* parent) { _parent=parent; }
			Node* getParent() const { return _parent; }
			virtual Tvaltype valueType() const =0;
			bool isInt() const { return valueType()==valtypeINT; };
			bool isFloat() const { return valueType()==valtypeFLOAT; };
			virtual Value evaluate() const =0;
			virtual Tfloat evaluateFloat() const { return Tfloat(evaluate()); }
			virtual Tint evaluateInt() const { return Tint(evaluate()); }
			virtual void write(std::ostream& ostrm) const =0;
			operator std::string() const
			{
				std::ostringstream oss;
				write(oss);
				return oss.str();
			}
			void write() const { write(std::cout); }
			void traverse(std::list<const Node*>& nodelist, enum Etraversetype traversetype=traversetypePOSTFIX) const;
			void writeSubExpr(std::ostream& ostrm=std::cout, enum Etraversetype traversetype=traversetypePOSTFIX, char sep=' ') const;

		protected:
			Node* _children[2];
			Node* _parent;
			short _priority;
	};
	// ------------------------------------------------------------------------------------ Node

	// NodeConstant ----------------------------------------------------------------------------
	class NodeConstant : public Node
	{
		public:
			NodeConstant(Tint val=0, Node* parent=NULL) : Node(NULL,NULL,parent,1), _value(val) {}
			NodeConstant(Tfloat val, Node* parent=NULL) : Node(NULL,NULL,parent,1), _value(val) {}
			Tvaltype valueType() const { return _value.getType(); }
			Value evaluate() const { return _value; }
			Tfloat evaluateFloat() const { return _value.getValueFloat(); }
			Tint evaluateInt() const { return _value.getValueInt(); }
			void write(std::ostream& ostrm) const { _value.write(ostrm); }
		protected:
			Value _value;
	};
	// ---------------------------------------------------------------------------- NodeConstant

	// NodeVariable ----------------------------------------------------------------------------
	class NodeVariable : public Node
	{
		public:
			NodeVariable(Variable* var, Node* parent=NULL) : Node(NULL,NULL,parent,1), _var(var) {}
			Tvaltype valueType() const { if(_var) return _var->getType(); else return valtypeNONE; }
			Value evaluate() const { if(_var) return _var->getValue(); else return Value(); }
			Tfloat evaluateFloat() const { if(_var) return _var->getValueFloat(); else return 0.; }
			Tint evaluateInt() const { if(_var) return _var->getValueInt(); else return 0; }
			void write(std::ostream& ostrm) const { if(_var) _var->write(ostrm); else ostrm << "undefVar"; }
		protected:
			Variable* _var;
	};
	// ---------------------------------------------------------------------------- NodeVariable

	// NodeOperation2 --------------------------------------------------------------------------
	class NodeOperation2 : public Node
	{
		public:
			NodeOperation2(char op, Node* child1, Node* child2, Node* parent=NULL);
			char op() const { return _operator; }
			Tvaltype valueType() const
			{
				if(_children[0] && _children[1])
				{
					if(_children[0]->isInt()&&_children[1]->isInt())
						return valtypeINT;
					else
						return valtypeFLOAT;
				}
				else
					return valtypeNONE;
			}
			Value evaluate() const;
			void write(std::ostream& ostrm) const { ostrm << _operator; }
		protected:
			char _operator;
	};
	// -------------------------------------------------------------------------- NodeOperation2

	// NodeFunction ----------------------------------------------------------------------------
	class NodeFunction : public Node
	{
		public:
			enum Efunctype {functypeNONE=0
			              , functypeMarker1Arg	// marker for functions with 1 argument ---
			              , functypeABS	// absolute value
			              , functypeFLOAT	// floating point value
			              , functypeINT	// integer value
			              , functypeFLOOR	// floor/round down integer value
			              , functypeCEIL	// ceiling/round up integer value
			              //, functypeTRUNC	// truncated integer value
			              , functypeROUND	// rounded integer value
			              , functypeSQRT	// square root
			              , functypeLN	// natural logarithm
			              , functypeLB	// binary logarithm
			              , functypeLG	// decimal logarithm
			              , functypeEXP	// exponential functions
			              , functypeSIN	// sine
			              , functypeCOS	// cosine
			              , functypeTAN	// tangent
			              , functypeASIN	// arc sine
			              , functypeACOS	// arc cosine
			              , functypeATAN	// arc tangent
			              , functypeMarker2Arg	// marker for functions with 2 arguments ---
			              , functypeMIN	// minimum of 2 values
			              , functypeMAX	// maximum of 2 values
			              , functypeMOD	// modulo function
			              , functypePOW	// power function
			              , functypeMarkerVarSet	// marker for functions using the VariableSet ===
			              , functypeMarkerVarSet1Arg	// marker for functions using the VariableSet with 1 argument ---
			              , functypeRCL	// recall stored value <id>
			              , functypeMarkerVarSet2Arg	// marker for functions using the VariableSet with 2 arguments ---
			              , functypeSTO	// store value to <id>
			               };

			static Efunctype functype(const std::string& name);

			NodeFunction(Efunctype func
			            ,Node* child1, Node* child0=NULL, Node* parent=NULL)
				: Node(child0,child1,parent,0), _functype(func) {}
			// functions with 1 argument use _children[1]
			// TODO: probably they should use _children[0] (again), but Node::writeSubExpr needs to be adapted
			//       also NodeFunctionStore
			Tvaltype valueType() const;
			Value evaluate() const;
			void write(std::ostream& ostrm) const;

		protected:
			enum Efunctype _functype;
	};
	// ---------------------------------------------------------------------------- NodeFunction

	// NodeFunctionVarSet ----------------------------------------------------------------------
	/*
		Functions with the capability to store and load values from the given VariableSet
	*/
	class NodeFunctionVarSet : public NodeFunction
	{
		public:
			// use enum Efunctype defined in NodeFunction

			NodeFunctionVarSet(Efunctype func, VariableSet *variableset
			                 ,Node* child1, Node* child0=NULL, Node* parent=NULL)
				: NodeFunction(func,child1,child0,parent), _variableset(variableset) {}
			// functions with 1 argument use _children[1] like NodeFunction (see comments there)
			Tvaltype valueType() const;
			Value evaluate() const;
			void write(std::ostream& ostrm) const;

		protected:
			VariableSet* _variableset;
	};
	// ---------------------------------------------------------------------- NodeFunctionVarSet

	// -------------------------------------------------------------------- Node and derivatives

		Expression(const std::string& label=std::string(), VariableSet* varset=NULL)
		          : _rootnode(NULL), _label(label), _variableset(varset), _variablesetcreated(false)
		{
			if(!_variableset)
			{
				_variableset=new VariableSet();
				_variablesetcreated=true;
			}
		}
		Expression(const Expression& expr) { *this=expr; }
		~Expression()
		{
			clear();
			if(_variablesetcreated) delete _variableset;
		}

		void clear() { if(_rootnode) { delete(_rootnode); _rootnode=NULL; } }
		void setLabel(const std::string& label) { _label=label; }
		const std::string& getLabel() const { return _label; }

		void initializeRPN(const std::string& exprstr, bool genlabel=true);

		Expression& operator=(const Expression& rhs )
		{
			_label=rhs._label;
			_variableset=rhs._variableset;
			_variablesetcreated=false;
			std::ostringstream oss;
			rhs.writeExpr(oss,traversetypePOSTFIX);
			_rootnode=NULL;	// prevents initializeRPN to clear a non-existant tree
			initializeRPN(oss.str());
			return *this;
		}

		bool isEmpty() const { return _rootnode==NULL; }
		bool isInt() const { if(_rootnode) return _rootnode->isInt(); else return false; }
		bool isFloat() const { if(_rootnode) return _rootnode->isFloat(); else return false; }
		Tfloat evaluateFloat() const
		{
			if(_rootnode)
				return _rootnode->evaluateFloat();
			else
				return 0.;
		}
		Tint evaluateInt() const
		{
			if(_rootnode)
				return _rootnode->evaluateInt();
			else
				return 0;
		}

		VariableSet* getVariableSet() { return _variableset; }
		Variable* getVariable(const std::string& name) { return _variableset->getVariable(name); }
		unsigned int VariablesCount() const { return _variableset->VariablesCount(); }
		bool existVariable(const std::string& name) const { return _variableset->existVariable(name); }
		unsigned int VariableGroupsCount() const { return _variableset->VariableGroupsCount(); }
		bool existVariableGroup(const std::string& name) const { return _variableset->existVariableGroup(name); }

		void writeExpr(std::ostream& ostrm=std::cout, enum Etraversetype traversetype=traversetypePOSTFIX, char sep=' ') const
		{
			if(_rootnode) _rootnode->writeSubExpr(ostrm,traversetype,sep);
		}
		operator std::string() const
		{
			std::ostringstream oss;
			writeExpr(oss,traversetypeINFIX,0);
			return oss.str();
		}
		void genLabel()
		{
			_label=operator std::string();
			//_label=static_cast<std::string>(*this);
		}

	protected:
		Node* _rootnode;
		std::string _label;
		VariableSet* _variableset;
		bool _variablesetcreated;
};
// ------------------------------------------------------------------------------ Expression


inline std::ostream& operator << (std::ostream& ostrm, const Expression::Value& v)
{
	v.write(ostrm);
	return ostrm;
}

inline std::ostream& operator << (std::ostream& ostrm, const Expression::Variable& v)
{
	v.write(ostrm);
	return ostrm;
}

inline std::ostream& operator << (std::ostream& ostrm, const Expression::Node& n)
{
	n.write(ostrm);
	return ostrm;
}

inline std::ostream& operator << (std::ostream& ostrm, const Expression& e)
{
	e.writeExpr(ostrm,Expression::traversetypeINFIX,0);
	return ostrm;
}


#endif
