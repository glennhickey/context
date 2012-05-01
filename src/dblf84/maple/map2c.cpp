#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cctype>

using namespace std;

struct Triple
{
	size_t x;
	size_t y;
	size_t z;
};


size_t scanBack(const char* buffer, size_t pos)
{
	int brack = 0;
	bool halt = false;
	while (!halt)
	{
		pos -= 1;
		switch(buffer[pos]) {
		case ')' : ++brack; break;
		case '(' : --brack;
		case '+' :
		case '-' :
		case '/' :
		case '*' : halt = brack <= 0; break;
		default : break;
		}
	}
	return pos + 1;
}

size_t scanFront(const char* buffer, size_t pos)
{
	int brack = 0;
	bool halt = false;
	while (!halt)
	{
		pos += 1;
		switch(buffer[pos]) {
		case '(' : ++brack; break;
		case ')' : --brack;
		case '+' :
		case '-' :
		case '/' :
		case '*' : halt = brack <= 0; break;
		default : break;
		}
	}
	return pos - 1;
}

void makeIndex(const char* buffer, vector<Triple>& carats)
{
	carats.clear();
	for (size_t i = 0; buffer[i] != '\0'; ++i)
	{
		if (buffer[i] == '^')
		{
			Triple trip;
			trip.x = scanBack(buffer, i);
			trip.y = i;
			trip.z = scanFront(buffer, i);
			carats.push_back(trip);

         if (trip.z == trip.y || trip.x == trip.y)
         {
            cerr << "BAD CARAT " << "pos=" << i  <<"  "
                 << trip.x << ", " << trip.y << ", " << trip.z << " " 
                 << string(buffer + trip.x, trip.z - trip.x + 1) << endl;
            exit(1);
         }
		}
	}
}

bool endOfFloat(const char* buffer, size_t pos)
{
   if (isdigit(buffer[pos]) && !isdigit(buffer[pos+1]) && buffer[pos+1] != ']')
   {
      while (isdigit(buffer[pos]) || buffer[pos] == 'e')
      {
         --pos;
      }
      if (buffer[pos] != '.')
      {
         return true;
      }
   }
   return false;
}

void convert(const char* buffer, const vector<Triple>& carats, ofstream& ofile)
{
	size_t cur = 0;
   size_t col = 15;
	for (size_t i = 0; buffer[i] != '\0'; ++i)
	{
		if (i >= carats[cur].x && i <= carats[cur].y)
		{
			string b(buffer + carats[cur].x, carats[cur].y - carats[cur].x);
			string t(buffer + carats[cur].y + 1, carats[cur].z - carats[cur].y);

         if (col + b.length() + t.length() + 6 > 80)
         {
            ofile << "\n";
            col = 0;
         }
   
         ofile << "pow(" << b << "," << t << ")";

         col += b.length() + t.length() + 6;
         
         if ((!isalnum(b[0]) && b[0] != '_') || (!isalnum(t[0]) && t[0] != '_'))
         {
            cerr << b << "  " << t << ",    ";
            cerr << string(buffer + carats[cur].x, carats[cur].z - carats[cur].x + 1) << endl;
         }
         
			i = carats[cur].z;
			++cur;
		}
		else
		{
         if (buffer[i] == '^')
         {
            cerr << "CARAT FOUND AT " << i << endl;
            cerr << "prev car " << carats[cur-1].x << "," << carats[cur-1].y << ","
                 << carats[cur-1].z << "  cur car " << carats[cur].x << "," << carats[cur].y
                 << "," << carats[cur].z << endl;
            cerr << string(buffer + carats[cur-1].x, carats[cur-1].z - carats[cur-1].x) << endl;
            exit(1);
         }
        
         if (col + 1 > 80 && (buffer[i] == '+' || buffer[i] == '-' || buffer[i] == '/'
                              || buffer[i] == '*'))
         {
            ofile << "\n";
            col = 0;
         }
         
			ofile << buffer[i];
         if (endOfFloat(buffer, i))
            ofile << ".";


         ++col;
		}
	}
}

int main(int argc, char** argv)
{
	if (argc != 3)
	{
		cerr << "map2c <input> <output>" <<endl;
		return 1;
	}

	size_t N = 100 * (1 << 20);

	char* buffer = new char[N];
   buffer[0] = '\0';

	ifstream ifile(argv[1]);
   while (ifile)
   {
      ifile.getline(buffer + strlen(buffer), N-1);
	}
	ofstream ofile(argv[2]);

	vector<Triple> carats;
	makeIndex(buffer, carats);
	
	ofile << "#include <cmath>\n"
			<< "#include \"dblhky.h\"\n\n"
			<< "using namespace std;\n\n"
			<< "void DblHKY::funname()\n{\n"
			<< "varname = ";
	
	convert(buffer, carats, ofile);

	ofile << ";\n}\n" <<endl;

	delete buffer;
	return 0;
}
