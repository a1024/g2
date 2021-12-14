//best viewed with tab size of 4 spaces
//g2_file.h - Include for file operations.
//Copyright (C) 2012-2020  Ayman Wagih Mohsen, unless source link provided.
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <https://www.gnu.org/licenses/>.
#ifndef			G2_FILE_H
#define			G2_FILE_H
#include		<vector>
#include		<string>
extern char		*programpath;//full path to program
extern const char *statedir,//full path to the folder containing the state, ends with back slash
				*statefilepath;//full path to file containing the state data
void			init_directories();

int				saveFile(const char *addr, const char *data, size_t size, bool binary=true);
int				loadFile(const char *addr, char *&data, size_t &size, bool binary=true);

void			directorycontents(const char *dirpath, std::vector<std::string> &ret);
void			mkdir_firsttime(const char *dirpath);
unsigned 		hexstr2uint(const char *text);
int				statefolder_readstate(char *&text, size_t &length);
void			statefolder_deletecontents();
const char*		version2str();
void 			statefolder_writestate(const char *usertext=nullptr);

const char*		loadresource(int name, int type, int &size);
#endif