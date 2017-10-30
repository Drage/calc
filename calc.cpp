
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stack>
#include <cmath>
#include <fstream>

using namespace std;

enum Associativity { LEFT, RIGHT, NONE };
const int MAX_ARGS = 2;
typedef double (*EvalFunc)(double*);
double x = 0;
typedef double (*EqFunc)(double, vector<string>, vector<string>);

struct Operator
{
	string value;
	int numArgs;
	EvalFunc evaluate;
	Associativity associativity;
	int precedence;

	Operator(string value, int numArgs, Associativity associativity, int precedence, EvalFunc func) 
	{ 
		this->value = value;
		this->numArgs = numArgs; 
		this->associativity = associativity;
		this->precedence = precedence;
		this->evaluate = func;
	}
};

vector<Operator*> operators;

Operator* GetOperator(string token)
{
	for (unsigned i = 0; i < operators.size(); i++)
	{
		if (operators[i]->value == token)
			return operators[i];
	}
	return NULL;
}

bool IsOperator(string token)
{
	return GetOperator(token) != NULL;
}

double Add(double* args) { return args[1] + args[0]; }
double Subtract(double* args) { return args[1] - args[0]; }
double Multiply(double* args) { return args[1] * args[0]; }
double Divide(double* args) { return args[1] / args[0]; }
double Power(double* args) { return pow(args[1], args[0]); }
double Sin(double* args) { return sin(args[0]); }
double Cos(double* args) { return cos(args[0]); }
double Tan(double* args) { return tan(args[0]); }
double Asin(double* args) { return asin(args[0]); }
double Acos(double* args) { return acos(args[0]); }
double Atan(double* args) { return atan(args[0]); }
double Sqrt(double* args) { return sqrt(args[0]); }
double Abs(double* args) { return abs(args[0]); }
double Log(double* args) { return log10(args[0]); }
double Ln(double* args) { return log(args[0]); }
double Pi(double* args) { return M_PI; }
double E(double* args) { return M_E; }
double X(double* args) { return x; }

double GetLastResult(double* args)
{
	// Open result file in temporary directory
	string ansFile = string(getenv("TMPDIR")) + "calc-result.bin";
	ifstream file(ansFile.c_str(), ios::in | ios::binary);
	if (file.is_open())
	{	
		// Read last result from file
		double lastResult;
		file.read((char*)&lastResult, sizeof(lastResult));
		file.close();
		return lastResult;
	}
	else
		throw string("previous result could not be loaded");
}

void SetLastResult(double result)
{
	// Open result file in temporary directory
	string resFile = string(getenv("TMPDIR")) + "calc-result.bin";
	ofstream file(resFile.c_str(), ios::out | ios::binary | ios::trunc);
	if (file.is_open())
	{
		// write result to file
		file.write((char*)&result, sizeof(result));
		file.close();
	}
	else
		throw string("could not save result");
}

void InitOperators()
{
	// Operators
	operators.push_back(new Operator("+", 2, LEFT, 2, Add));
	operators.push_back(new Operator("-", 2, LEFT, 2, Subtract));
	operators.push_back(new Operator("*", 2, LEFT, 3, Multiply));
	operators.push_back(new Operator("/", 2, LEFT, 3, Divide));
	operators.push_back(new Operator("^", 2, RIGHT, 4, Power));
	
	// Functions
	operators.push_back(new Operator("sin", 1, NONE, 5, Sin));
	operators.push_back(new Operator("cos", 1, NONE, 5, Cos));
	operators.push_back(new Operator("tan", 1, NONE, 5, Tan));
	operators.push_back(new Operator("asin", 1, NONE, 5, Asin));
	operators.push_back(new Operator("acos", 1, NONE, 5, Acos));
	operators.push_back(new Operator("atan", 1, NONE, 5, Atan));
	operators.push_back(new Operator("sqrt", 1, NONE, 5, Sqrt));
	operators.push_back(new Operator("abs", 1, NONE, 5, Abs));
	operators.push_back(new Operator("log", 1, NONE, 5, Log));
	operators.push_back(new Operator("ln", 1, NONE, 5, Ln));
	
	// Constants
	operators.push_back(new Operator("pi", 0, NONE, 5, Pi));
	operators.push_back(new Operator("e", 0, NONE, 5, E));
	operators.push_back(new Operator("ans", 0, NONE, 5, GetLastResult));
	operators.push_back(new Operator("x", 0, NONE, 5, X));
	
	// Brackets
	operators.push_back(new Operator("(", 0, NONE, 0, NULL));
	operators.push_back(new Operator(")", 0, NONE, 0, NULL));
}

// Parse expression into list of operators and operands
void Tokenise(string expression, vector<string> &tokens)
{
	while (expression.length() > 0)
	{
		if (isdigit(expression[0]))
		{
			// Find the string index where the number ends
			int numEnd = 1;
			while (numEnd < expression.length() && (isdigit(expression[numEnd]) || expression[numEnd] == '.'))
				numEnd++;
			
			// Add number to token list and remove from expression
			tokens.push_back(expression.substr(0, numEnd));
			expression = expression.substr(numEnd);
		}
		else if (isalpha(expression[0]))
		{
			// Find the string index where the function ends
			int funcEnd = 1;
			while (funcEnd < expression.length() && isalpha(expression[funcEnd]))
				funcEnd++;
			
			// Check function is a legit operator (sin, cos, tan, sqrt or pi)
			if (IsOperator(expression.substr(0, funcEnd)))
			{
				// Add function to token list and remove from expression
				tokens.push_back(expression.substr(0, funcEnd));
				expression = expression.substr(funcEnd);
			}
			else
			{
				stringstream ss;
				ss << "invalid function '" << expression.substr(0, funcEnd) << "'";
				throw ss.str();
			}
		}
		else if (IsOperator(expression.substr(0, 1)))
		{
			// Add operator to token list and remove from expression
			// (All operators which are not functions are 1 character)
			tokens.push_back(expression.substr(0, 1));
			expression = expression.substr(1);
		}
		else if (expression[0] == ' ')
		{
			// Ignore spaces
			expression = expression.substr(1);
		}
		else
		{
			// Unknown character
			stringstream ss;
			ss << "invalid character '" << expression[0] << "'";
			throw ss.str();
		}
	}
	
	// Handle negative numbers
	// For any sequence [operator, '-', number] make the number -ve and remove the '-' operator
	for (int i = 0; i < (int)tokens.size() - 2; i++)
	{
		if ((tokens[i] == "+" || tokens[i] == "-" || tokens[i] == "*" || tokens[i] == "/")
			&& tokens[i+1] == "-" && isdigit(tokens[i+2][0]))
		{
			tokens[i+2] = "-" + tokens[i+2];
			tokens.erase(tokens.begin()+i+1);
		}
	}
	// Special case where first token is '-' operator
	if (tokens[0] == "-")
	{
		if (isdigit(tokens[1][0]))
		{
			tokens[1] = "-" + tokens[1];
			tokens.erase(tokens.begin());
		}
	}
	
	// Handle implicit multiplication
	// E.g. 2(4-1) should be changed to 2*(4-1)
	for (int i = 0; i < (int)tokens.size() - 1; i++)
	{
		if (isdigit(*tokens[i].rbegin()) && (tokens[i+1] == "(" || tokens[i+1] == "x"))
			tokens.insert(tokens.begin()+i+1, "*");
	}
	// Special case where first token is ')' and second is 'x'
	// i.e. (1/2)x
	for (int i = 0; i < (int)tokens.size() - 1; i++)
	{
		if (tokens[i] == ")" && tokens[i+1] == "x")
			tokens.insert(tokens.begin()+i+1, "*");
	}
}

// Convert infix notation to reverse polish notation using the shunting-yard algorithm
void ConvertToRPN(vector<string> &infix, vector<string> &rpn)
{
	stack<Operator*> opStack;

	for (vector<string>::iterator i = infix.begin(); i != infix.end(); i++)
	{
		Operator *op = GetOperator(*i);
		
		if (op == NULL)
		{
			// Number
			rpn.push_back(*i);
		}
		else
		{
			// Brackets
			if (op->value == "(")
				opStack.push(op);
			else if (op->value == ")")
			{
				while (opStack.top()->value != "(")
				{
					rpn.push_back(opStack.top()->value);
					opStack.pop();

					if (opStack.size() == 0)
						throw "bracket mismatch";
				}
				opStack.pop();
			}
			else
			{
				// Operator
				while (opStack.size() > 0)
				{
					if ((op->associativity == LEFT && op->precedence <= opStack.top()->precedence) 
						|| op->precedence < opStack.top()->precedence)
					{
						rpn.push_back(opStack.top()->value);
						opStack.pop();
					}
					else
						break;
				}
				opStack.push(op);
			}
		}
	}
	while (opStack.size() != 0)
	{
		if (opStack.top()->value == "(")
			throw "bracket mismatch";

		rpn.push_back(opStack.top()->value);
		opStack.pop();
	}
}

// Evaluate an expression in RPN form
double Evaluate(vector<string> &tokens)
{
	stack<double> rpnStack;
	for (unsigned i = 0; i < tokens.size(); i++)
	{
		if (!IsOperator(tokens[i]))
		{
			// Convert to double and push onto stack
			double value;
			stringstream converter(tokens[i]);
			converter >> value;
			rpnStack.push(value);
		}
		else
		{
			Operator *op = GetOperator(tokens[i]);
			
			// Check that there are sufficient operands on the stack
			if ((int)rpnStack.size() < op->numArgs)
				throw string("failed to evaluate expression");
			
			// Create args list
			double args[MAX_ARGS];
			for (int n = 0; n < op->numArgs; n++)
			{
				args[n] = rpnStack.top();
				rpnStack.pop();
			}
			
			// Evaluate operator
			double result = op->evaluate(args);
			
			// Put result back on stack
			rpnStack.push(result);
		}
	}
	
	// Check only one number is left on the stack
	if (rpnStack.size() != 1)
		throw string("failed to evaluate expression");
	
	return rpnStack.top();
}

// Implementation of Brent's algorithm for root-finding
// (a combination of bisection, secant, and inverse quadratic interpolation)
double Solve(EqFunc function, vector<string> lhs, vector<string> rhs, double lowerLimit, double upperLimit, double errorTolerance)
{
	double a = lowerLimit;
	double b = upperLimit;
	double c = 0;
	double d = numeric_limits<double>::max();
	
	double fa = function(a, lhs, rhs);
	double fb = function(b, lhs, rhs);
	
	double fc = 0;
	double s = 0;
	double fs = 0;
	
	// Swap a and b if |f(a)| < |f(b)|
	if (abs(fa) < abs(fb))
	{ 
		double tmp = a; 
		a = b; 
		b = tmp; 
		tmp = fa; 
		fa = fb; 
		fb = tmp; 
	}
	
	c = a;
	fc = fa;
	bool mflag = true;
	int i = 0;
	
	while (fb != 0 && abs(a-b) > errorTolerance)
	{
		if (fa != fc && fb != fc)
		{
			// Inverse quadratic interpolation
			s = a * fb * fc / (fa - fb) / (fa - fc) 
				+ b * fa * fc / (fb - fa) / (fb - fc) 
				+ c * fa * fb / (fc - fa) / (fc - fb);
		}
		else
		{
			// Secant Rule
			s = b - fb * (b - a) / (fb - fa);
		}
		
		// Bisection method
		double tmp = (3 * a + b) / 4;
		if (!((s > tmp && s < b) 
			|| (s < tmp && s > b)) 
			|| (mflag && abs(s - b) >= abs(b - c) / 2) 
			|| (!mflag && abs(s - b) >= abs(c - d) / 2))
		{
			s = (a + b) / 2;
			mflag = true;
		}
		else
			mflag = false;
		
		fs = function(s, lhs, rhs);
		d = c;
		c = b;
		fc = fb;
		if (fa * fs < 0) 
		{ 
			b = s; 
			fb = fs; 
		}
		else 
		{ 
			a = s; 
			fa = fs; 
		}
		
		// Swap (a,b) if |f(a)| < |f(b)|
		if (abs(fa) < abs(fb))
		{ 
			double tmp = a; 
			a = b; 
			b = tmp; 
			tmp = fa; 
			fa = fb; 
			fb = tmp; 
		}
		
		// Keep track of iteration count, max out at 1000
		i++;
		if (i > 1000)
			throw string("failed to solve equation");
	}
	
	// return root
	return b;
}

// Evaluate the equation for a given value
double SolveTest(double value, vector<string> lhs, vector<string> rhs)
{
	x = value;
	return Evaluate(lhs) - Evaluate(rhs);
}

void ProcessInput(string input) {
	try
	{
		// Check if simple expression or equation
		int equals = input.find("=");
		if (equals == string::npos) 
		{
			// Break up into tokens
			vector<string> infix;
			Tokenise(input, infix);
			
			// Convert to RPN
			vector<string> tokens;
			ConvertToRPN(infix, tokens);
		
			// Evaluate expression
			double result = Evaluate(tokens);
			
			// Save result so it can be referenced by subsequent calculations
			SetLastResult(result);
			
			// Display result
			cout << "= " << result << endl;
		}
		else
		{
			// Split string into lhs and rhs of equation
			string inputL = input.substr(0, equals);
			string inputR = input.substr(equals + 1, input.length()-equals-1);

			// Break up into tokens
			vector<string> infixL, infixR;
			Tokenise(inputL, infixL);
			Tokenise(inputR, infixR);

			// Convert to RPN
			vector<string> lhs, rhs;
			ConvertToRPN(infixL, lhs);
			ConvertToRPN(infixR, rhs);

			// Solve equation
			double x = Solve(SolveTest, lhs, rhs, -1000, 1000, 0.00001);

			// Save result so it can be referenced by subsequent calculations
			SetLastResult(x);

			// Display result
			cout << "x = " << x << endl;
		}
	}
	catch (string errorMsg) 
	{
		cout << "Error: " << errorMsg << endl;
	}
}

int main(int argc, char **argv)
{
	// Create operator list
	InitOperators();

	// Get expression from arguments
	if (argc == 2)
		ProcessInput(argv[1]);
	else
	{
		for (;;)
		{
			// Prompt user to enter expression
			cout << "> ";
			char buffer[512];
			cin.getline(buffer, sizeof(buffer));
			string expression = buffer;
	
			// End input loop
			if (expression == "exit")
				break;
	
			ProcessInput(expression);
		}
	}

	return 0;
}
