#pragma once

#include "pdbchain.h"
#include "seqdb.h"
#include "arrays.h"

class CMProf
	{
public:
	SquareMatrix<double> m_MeanDistMx;
	SquareMatrix<double> m_StdDevs;
	const SeqDB *m_MSA = 0;

// Training only
	vector<bool> m_ColIsCore;
	vector<uint> m_CoreCols;
	vector<SquareMatrix<double> > m_DistMxVec;
	map<string, uint> m_UngappedSeqToIdx;


public:
	CMProf()
		{
		Clear();
		};

	void Clear()
		{
		m_ColIsCore.clear();
		m_MeanDistMx = SquareMatrix<double>();
		m_StdDevs = SquareMatrix<double>();
		m_DistMxVec.clear();
		}
	
	void SetMSA(const SeqDB &MSA);
	uint GetColCount() const { return m_MSA->GetColCount(); }
	uint GetCoreColCount() const { return SIZE(m_CoreCols); }

	void ToFile(const string &FileName) const;
	void MxToFile(FILE *f, const string &Name,
	  const SquareMatrix<double> &Mx) const;
	void MxFromFile(FILE *f, string &Name, uint CoreColCount,
	  SquareMatrix<double> &Mx);
	void FromFile(FILE *f);
	void FromFile(const string &FileName);

// For training
public:
	void InitTrain()
		{
		m_DistMxVec.clear();
		}

	void FinalizeTrain();
	bool TrainChain(const PDBChain &Chain);
	void GetDistMx(const PDBChain &Chain, const vector<uint> &PosVec,
	  SquareMatrix<double> &DistMx);
	void GetMeanStdDev(uint i, uint j,
	  double &Mean, double &StdDev) const;
	};

double GetNormal(double Mu, double Sigma, double x);
