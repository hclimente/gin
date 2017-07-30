//
// Created by hclimente on 25/07/2017.
//

#ifndef GIN_SETTINGS_H
#define GIN_SETTINGS_H

class Settings
{
public:
	static S& getInstance()
	{
		static S    instance; // Guaranteed to be destroyed.
		// Instantiated on first use.
		return instance;
	}
private:
	S() {}                    // Constructor? (the {} brackets) are needed here.

	string __pedBasename;
	string __pedBasename;

public:
	S(S const&)               = delete;
	void operator=(S const&)  = delete;

};

#endif //GIN_SETTINGS_H
