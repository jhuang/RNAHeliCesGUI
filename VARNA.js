/**
 * @author ponty
 */
function setTitle(appletid,ntitle)
{
    var applet = document.getElementById(appletid);
    var script = "setTitle(\""+ntitle+"\")";
	applet.runScript(script);
};

function setSeq(appletid,nseq)
{
    var applet = document.getElementById(appletid);
    var script = "setSeq(\""+nseq+"\")";
	applet.runScript(script);
};

function setStruct(appletid,nstr)
{
    var applet = document.getElementById(appletid);
    var script = "setStruct(\""+nstr+"\")";
	applet.runScript(script);
};

function setStructSmooth(appletid,nstr)
{
    var applet = document.getElementById(appletid);
    var script = "setStructSmooth(\""+nstr+"\")";
	applet.runScript(script);
};

function setRNA(appletid,nseq,nstr)
{
    var applet = document.getElementById(appletid);
    var script = "setRNA(\""+nseq+"\",\""+nstr+"\")";
	applet.runScript(script);
};
		

function setRNASmooth(appletid,nseq,nstr)
{
    var applet = document.getElementById(appletid);
    var script = "setRNASmooth(\""+nseq+"\",\""+nstr+"\")";
	applet.runScript(script);
};
		
function redraw(appletid,nalgo)
{
    var applet = document.getElementById(appletid);
    var script = "redraw(\""+nalgo+"\")";
	applet.runScript(script);
};

function setColorMapValues(appletid,values)
{
    var applet = document.getElementById(appletid);
	var txt = "";
	for(var i=0;i<values.length;i++)
	{
		if (i>0)
		  txt += ", ";
		txt += values[i];
	}
    var script = "setValues(["+txt+"])";
	applet.runScript(script);
};
