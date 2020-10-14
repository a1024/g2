//best viewed with tab size of 4 spaces
//g2_file.cpp - Implementation of file operations.
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
#include		"g2_file.h"
#include		"g2_error.h"
#include		"g2_expr.h"
#include		<Windows.h>
#include		<sys/stat.h>
#include		<vector>
#include		<string>
char			*programpath=nullptr;//full path to program		//TODO: UNICODE
const char		statefoldername[]="g2_state",
				statefilename[]="state.txt",
				*statedir=nullptr,//full path to the folder containing the state, ends with back slash
				*statefilepath=nullptr;//full path to file containing the state data
const char*		alloc_str(std::string const &str)
{
	unsigned size=str.size();
	char *ret=(char*)malloc(size+1);
	memcpy(ret, str.c_str(), size);
	ret[size]='\0';
	return ret;
}
#define 		ALLOC_STR(pointer, str)		if(!(pointer))pointer=alloc_str(str)
void			init_directories()
{
	if(!programpath)
	{
		programpath=(char*)malloc(MAX_PATH+1);
		GetModuleFileNameA(nullptr, programpath, MAX_PATH);
		for(int k=strlen(programpath)-1;k>=0;--k)//remove executable name
		{
			if(programpath[k]=='/'||programpath[k]=='\\')
			{
				programpath[k+1]='\0';
				break;
			}
		}
	}
	ALLOC_STR(statedir, std::string(programpath)+statefoldername+'\\');
	ALLOC_STR(statefilepath, std::string(statedir)+statefilename);
}

int				saveFile(const char *addr, const char *data, size_t size, bool binary)
{
	const char mode[]={'w', char(binary*'b'), '\0'};
	FILE *file;
	int error=fopen_s(&file, addr, mode);	//SYS_CHECK();//2: No such file or directory
//	FILE *file=fopen(addr, binary?"wb":"w");	SYS_CHECK();
	if(error)
		LOGERROR("Failed to save %s", addr);
	else
	{
		int byteswritten=fwrite(data, 1, size, file);	//SYS_CHECK();//2: No such file or directory
		if(byteswritten!=size)
			LOGERROR("saved %d/%d bytes of %s", byteswritten, size, addr);
		//	LOGERROR("Failed to save %s", addr);
		int error=fclose(file);
		if(error==EOF)
			SYS_CHECK();
		return 0;
	}
#if 0
	std::ofstream file(addr, std::ios::out|std::ios::binary);
	if(file.is_open())
	{
		file.write(data, size);
		file.close();
		return 0;
	}
#endif
	return 1;
}
long			getfilesize(FILE *file)
{//https://www.cplusplus.com/reference/cstdio/fread/
	fseek(file , 0 , SEEK_END);	SYS_CHECK();
	long size=ftell(file);
	rewind(file);
	return size;
}
int				loadFile(const char *addr, char *&data, size_t &size, bool binary)
{
#if 0
	FILE *file=fopen(addr, binary?"rb":"r");
	int error=errno;
	const char *msg=strerror(error);
	if(file)
	{
		size=getfilesize(file);
		data=(char*)malloc(size);
		size_t bytesread=fread(data, 1, size, file);	SYS_CHECK();
		if(bytesread!=size)
			LOGERROR("Error reading %s", addr);
		int error=fclose(file);	SYS_CHECK();
		if(error==EOF)
			SYS_CHECK();
		return 0;
	}
	size=0, data=nullptr;
	return 1;
#endif
#if 1
	struct stat info={};
	if(!stat(addr, &info))
	{
		SYS_CHECK();
		size=info.st_size;
		data=(char*)malloc((size+1)*sizeof(unsigned char));
		data[size]='\0';
		const char mode[]={'r', char(binary*'b'), '\0'};
		FILE *file;
		int error=fopen_s(&file, addr, mode);
	//	FILE *file=fopen(addr, binary?"rb":"r");
		if(error)
		{
			SYS_CHECK();
			free(data);
		}
		else
		{
			size_t bytesread=fread(data, 1, size, file);	SYS_CHECK();
			data[bytesread]='\0';
			int nnl=0;//number of newlines
			if(!binary)
				for(unsigned k=0;k<bytesread;++k)
					nnl+=data[k]=='\n';
			if(bytesread+nnl!=size)//Windows newlines "\r\n" (2 chars) read as "\n" (1 char)
				LOGERROR("Read %d/%d from %s", (int)bytesread, (int)size, addr);
			int error=fclose(file);
			if(error==EOF)
				SYS_CHECK();
			return 0;
		}
		//std::ifstream file(addr, std::ios::in|std::ios::binary);
		//if(file.is_open())
		//{
		//	file.read((char*)data, size);
		//	file.close();
		//	return 0;
		//}
	}
	size=0, data=nullptr;
	return 1;
#endif
}

void			directorycontents(const char *dirpath, std::vector<std::string> &ret)
{
	_WIN32_FIND_DATAA data;
	std::string dir=dirpath;
	dir+='*';
	//char &c=*dir.rbegin();
	//if(c=='\\'||c=='/')
	//	dir.pop_back();
	void *hSearch=FindFirstFileA(dir.c_str(), &data);
	if(hSearch!=INVALID_HANDLE_VALUE)
	{
		do
		{
			if(strcmp(data.cFileName, ".")&&strcmp(data.cFileName, ".."))
				ret.push_back(data.cFileName);
		}
		while(FindNextFileA(hSearch, &data));
		FindClose(hSearch);
	}
}
void			mkdir_firsttime(const char *dirpath)
{
	int attrib=GetFileAttributesA(dirpath);//https://stackoverflow.com/questions/3828835/how-can-we-check-if-a-file-exists-or-not-using-win32-program
	//struct stat s={};
	//int doesntexist=stat(dirpath, &s);//returns -1 even if folder exists
	if(attrib==INVALID_FILE_ATTRIBUTES)//
	{
		int success=CreateDirectoryA(dirpath, nullptr);
	//	auto fail=mkdir(dirpath, 0777);
		if(!success)
			SYS_CHECK();
	}
	else if(!(attrib&FILE_ATTRIBUTE_DIRECTORY))//a FILE exists with same name
		LOGERROR("Please put the program away from the \'g2_state\' file.");
}
unsigned 		hexstr2uint(const char *text)
{
	unsigned number=0;
	if(text)
	{
		for(unsigned k=0;k<8;++k)
		{
			unsigned nibble=0;
			if(text[k]>='0'&&text[k]<='9')
				nibble=char(text[k]-'0');
			else if(text[k]>='A'&&text[k]<='F')
				nibble=char(10+text[k]-'A');
			else if(text[k]>='a'&&text[k]<='f')
				nibble=char(10+text[k]-'a');
			else
				break;
			number|=nibble<<((7-k)<<2u);
		}
	}
	return number;
}
int				statefolder_readstate(char *&text, size_t &length)
{
	init_directories();
	mkdir_firsttime(statedir);
	return loadFile(statefilepath, text, length, false);//zero if found (success)
}
void			statefolder_deletecontents()
{
	init_directories();
	std::vector<std::string> contents;
	directorycontents(statedir, contents);
	for(int k=0, nitems=contents.size();k<nitems;++k)
	{
		const char *name=contents[k].c_str();
	//	if(strcmp(name, ".")&&strcmp(name, ".."))
	//	{
		//	LOGI("Deleting %s", name);
			int error=remove((statedir+('/'+contents[k])).c_str());
			if(error)
				SYS_CHECK();
	//	}
	}
}
const char*		version2str()
{
	static char vstr[9]={};
	for(unsigned k=0;k<8;++k)
	{
		char nibble=char(g2_version.id>>((7-k)<<2u)&15u);
		vstr[k]=char((nibble>=10?'A':'0')+nibble%10);
	}
	return vstr;
}
void 			statefolder_writestate(const char *usertext)
{
	std::string str=version2str();
	str+='\n';
	if(usertext)
		str+=usertext;

	init_directories();
	saveFile(statefilepath, str.c_str(), str.size(), false);
}