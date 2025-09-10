#include "myutils.h"
#include "dssaligner.h"
#include "arrays.h"
#include "chainreader2.h"

void cmd_alignselfrev()
	{
	const string &QFN = g_Arg1;
	FILE *fOut = CreateStdioFile(opt(output));

	DSSParams Params;
	Params.SetDSSParams(DM_AlwaysSensitive);
	Params.m_UsePara = false;
	Params.m_Omega = 0;

	DSSAligner DA;
	DA.SetParams(Params);

	DSS D;
	D.SetParams(Params);

	ChainReader2 CR;
	CR.Open(QFN);
	uint N = 0;
	for (;;)
		{
		PDBChain *Chain = CR.GetNext();
		if (Chain == 0)
			break;
		++N;
		if (N%100 == 0)
			Progress("%u\r", N);
		PDBChain Rev = *Chain;
		Rev.Reverse();

		D.Init(*Chain);
		Matrix<byte> Profile = D.GetProfile();

		D.Init(Rev);
		Matrix<byte> RevProfile = D.GetProfile();

		DA.SetQuery(*Chain, &Profile, 0, 0, FLT_MAX);
		DA.SetTarget(Rev, &RevProfile, 0, 0, FLT_MAX);
		DA.AlignQueryTarget();
		DA.ToTsv(fOut, true);
		}
	CloseStdioFile(fOut);
	}
