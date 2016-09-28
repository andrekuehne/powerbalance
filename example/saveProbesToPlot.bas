'#Language "WWB-COM"

Option Explicit

Sub Main

' init variables
Dim CurrProbeTreeItem As String
Dim CurrProbeFileName As String
Dim FinalResultVoltages As Object
Dim FinalResultCurrents As Object
Dim FinalSMat As Object
Dim CurrProbeResult As Object
Dim NumOfResPoints As Long
Dim CurrProbeRealVal As Double
Dim CurrProbeImagVal As Double
Dim CurrProbeNum As Long
Dim CurrTxChNum As Long

' init final result
Set FinalResultVoltages = Result1DComplex("")
Set FinalResultCurrents = Result1DComplex("")
Set FinalSMat = Result1DComplex("")
Set CurrProbeResult = DS.Result1DComplex("")

' dialog for number of Tx channels
	Begin Dialog UserDialog 400,105,"Enter number of Tx channels" ' %GRID:10,7,1,1
		OKButton 60,70,90,21
		CancelButton 230,70,90,21
		TextBox 150,28,90,21,.numTxCh
	End Dialog
	Dim dlgTxCh As UserDialog
	dlgTxCh.numTxCh = "2"
	Dialog dlgTxCh

	Begin Dialog UserDialog1 400,105,"Enter number of probes" ' %GRID:10,7,1,1
		OKButton 60,70,90,21
		CancelButton 230,70,90,21
		TextBox 150,28,90,21,.numLumped
	End Dialog
	Dim dlgLumped As UserDialog1
	dlgLumped.numLumped = "9"
	Dialog dlgLumped



For CurrTxChNum = 1 To Val(dlgTxCh.numTxCh) STEP 1
	For CurrProbeNum = 1 To Val(dlgLumped.numLumped) STEP 1

		'Debug.Print CurrProbeNum

		' get voltage probe filename and create object
		CurrProbeTreeItem = "Tasks\SS_TxCh" & CStr(CurrTxChNum) & "\FD Voltages\P" & CStr(CurrProbeNum)
		CurrProbeFileName = DSResultTree.GetFileFromTreeItem(CurrProbeTreeItem)

		' load data into current probe result
		CurrProbeResult.Load(CurrProbeFileName)
		NumOfResPoints = CurrProbeResult.GetN

		' extract real and imaginbary parts
		CurrProbeRealVal = CurrProbeResult.GetYRe(0)
		CurrProbeImagVal = CurrProbeResult.GetYIm(0)

		FinalResultVoltages.AppendXY(CurrTxChNum,CurrProbeRealVal,CurrProbeImagVal)

		' get current probe filename and create object
		CurrProbeTreeItem = "Tasks\SS_TxCh" & CStr(CurrTxChNum) & "\FD Currents\P" & CStr(CurrProbeNum)
		CurrProbeFileName = DSResultTree.GetFileFromTreeItem(CurrProbeTreeItem)

		' load data into current probe result
		CurrProbeResult.Load(CurrProbeFileName)
		NumOfResPoints = CurrProbeResult.GetN

		' extract real and imaginbary parts
		CurrProbeRealVal = CurrProbeResult.GetYRe(0)
		CurrProbeImagVal = CurrProbeResult.GetYIm(0)

		FinalResultCurrents.AppendXY(CurrTxChNum,CurrProbeRealVal,CurrProbeImagVal)

	Next CurrProbeNum
Next CurrTxChNum


' Add to result
With FinalResultVoltages

    .Title("Lumped_Element_Voltages")
    .save("lumped_voltages.sig")
    .AddToTree("DS_Results\Lumped_Elements\Voltages")

End With

With FinalResultCurrents


    .Title("Lumped_Element_Currents")
    .save("lumped_currents.sig")
    .AddToTree("DS_Results\Lumped_Elements\Currents")

End With



End Sub
