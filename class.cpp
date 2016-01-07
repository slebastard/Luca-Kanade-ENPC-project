#include "options.h"
using namespace command;

exception no_value("Aucune valeur n'a �t� correctement renseign�e pour ce param�tre"), wrong_syntax("L'argument sp�cifi� est syntaxiquement incorrect"), unknown("L'argument pass� ne correspond � aucune option connue");


/* Methodes de OptionsField */

template <typename T>
OptionField<T>::OptionField(char *name[])
{
	opt_name = name;
}

template <typename T>
void OptionField<T>::filter(int& argc, char *argv[])
{
	if (argc > 1)
	{
		for (int rank = 1; rank < argc; rank++)
		{
			if (argv[rank][0] == '-' && argv[rank][1] == opt_name[0][0])	 // Attention: ne marche que si les options/arguments contiennent au moins 2 caract�res
			{
				if (rank + 1 < argc && argv[rank + 1][0] != '-')
					*arg = argv[rank + 1];
				else
					throw(no_value);
			}
			else if (strtok(argv[rank], '=') == opt_name[1])
			{
				if (char2string(last_part(argv[rank], '=')).size() != 0)
					*arg = last_part(argv[rank], '=');
				else
					throw(no_value);

			}
		}
	}
}

/* Sp�cialisation de OptionField poru des bool�ens */

template <>
OptionField<bool>::OptionField(char *name[])
{
	opt_name[0] = name[0];
	opt_name[1] = name[1];
	*arg = false;
}

template <>
void OptionField<bool>::filter(int& argc, char *argv[])
{
	if (argc > 1)
	{
		for (int rank = 1; rank < argc; rank++)
		{
			if (argv[rank][0] == '-' && argv[rank][1] != '-')
			{
				string str_arg = char2string(argv[rank]);
				if (str_arg.find_first_of(opt_name[0]))
					*arg = true;
			}
			else if (argv[rank] == opt_name[1])
				*arg = true;
		}
	}
}

/* M�thodes de CmdLine */
CmdLine::CmdLine() : opt()
{
}

CmdLine::CmdLine(vector<Option*> opt0) : opt()
{
	opt = opt0;
}

void CmdLine::add(Option *opt0)
{
	opt.push_back(opt0);
}

void CmdLine::process(int& argc, char *argv[])
{
		for (int rk = 1; rk < argc; rk++)
		{
			if (argv[rk][0] != '-' || (char2string(argv[rk]).size() == 1 && argv[rk][0] == '-')) // Erreur de syntaxe sur l'option
				throw wrong_syntax;
			else if (char2string(argv[rk]).size() > 1 && argv[rk][0] == '-') // Les param�tres commenceant par -, avec autre chose derri�re
			{
				/*if (char2string(argv[rk]).size() > 2 && argv[rk][1] != '-') // On v�rifie si le param�tre renseign� existe
				{
					bool exists = false;
					for (vector<Option*>::iterator opt_i = opt.begin(); opt_i != opt.end(); opt_i++)
						exists |= (argv[rk][1] == (*opt_i)->opt_name[0][0]);
					if (!exists)
						throw unknown;
				} */
			}
		}
	for (vector<Option*>::iterator opt_i = opt.begin(); opt_i != opt.end(); opt_i++)
		(*opt_i)->filter(argc, argv);
}

/* Autres fonctions utilises ici */

template<typename T>
Option* command::make_option(char c, string s1, T var)
{
	char* C = new char[1];
	*C = c;
	char *name[] = new char*[2];
	name[0] = C;
	name[1] = ("--" + s1).c_str();
	OptionField<T> *option = new OptionField<T>(name);
	return option;
}

string char2string(char* c)
{
	stringstream ss;
	string s;
	ss << c;
	ss >> s;
	return s;
}

string last_part(char *c, char limit)
{
	/* A compl�ter: doit retourner la partie droite (ie l'argument) dans
	     --year=2019
	*/
	string s = char2string(c);
	char* parsed;
	int limit_rk = s.size() - 1;
	for (int rk = s.size() - 1; rk != -1; rk--)
	{
		if (rk == limit)
		{
			limit_rk = rk;
			break;
		}
	}
	parsed = new char[s.size() - 1 - limit_rk];
	for (int rk = 0; rk != s.size() - 1 - limit_rk; rk++)
		parsed[rk] = c[limit_rk + 1 + rk];
	return char2string(parsed);
}