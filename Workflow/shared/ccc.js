// JScript source code

// parameter values used
var strTotalDays = "TotalDays";
var strUsedDays = "NumDays";

// rgb colour values
var strRedRGB = "rgb(184,29,3)"
var strGrnRGB = "rgb(31,156,35)";
var strOrangeRGB = "rgb(210,98,0)";

function reLayout(strURL)
{
	var nTrialTotal = 0;
	var bShowCounter = false;
	
	var strShow = extractParam(strURL, "ShowCounter");
	if(strShow.toLowerCase() == "yes")
	{
	    bShowCounter = true;
	}
	if ( bShowCounter )
	{
	    // get the values for the trial remaining
		var strTrialTotal = extractParam(strURL,strTotalDays);
		
		// get the values for the trial total
		if(strTrialTotal.length  > 0)
		{
			nTrialTotal = parseInt(strTrialTotal,10);
		}
		if ( nTrialTotal <= 0)
		{
		    bShowCounter = false;
		}
	}
	
	if ( bShowCounter)
	{
	    // only do this if the counter is to be shown   
		var tdPromo = document.getElementById("tdPromo")
		if (tdPromo != null)
		{
			tdPromo.style.display = "";
		}
		var divTrialRemain = document.getElementById("divTrialRemainCounter")
		if (divTrialRemain != null)
		{
			divTrialRemain.style.display = "";
		}
		var tdProgress = document.getElementById("tdProgress")
		if (tdProgress != null)
		{
			tdProgress.style.display = "";
		}	
			
		// get the values for the trial remaining
		var strDaysRemain = extractParam(strURL,strUsedDays);
		
		var nDaysRemain = 0;
		if(strDaysRemain.length > 0)
		{
			nDaysRemain = parseInt(strDaysRemain,10);
		}
		if (nDaysRemain < 0)
		{
			nDaysRemain = 0;
		}
		
		// insert the days remaining into the document
		var elmCounter = document.getElementById("divCountOnBanner");
		if(null != elmCounter)
		{
			elmCounter.innerText = String(nDaysRemain);
		}
		
		// adust the size of the progress bar
		var elmTopSlider = document.getElementById("divOverSlider");
		var elmBotSlider = document.getElementById("divUnderSlider");	
		var elmPImgLeft = document.getElementById("divPBStart");
		var elmPImgRight = document.getElementById("divPBEnd");
		// we now have to change the colour of several elements as well as swap some graphics
		var strColourRGB = strGrnRGB;	// initial RGB colour
		
		if( nDaysRemain <= 5 )
		{
			strColourRGB = strRedRGB;
			// red progressbar
			elmTopSlider.style.backgroundImage = "url(./shared/Images/TrialCounter/pbred.gif)";
			elmPImgLeft.style.backgroundImage = "url(./shared/Images/TrialCounter/pb_startred.gif)";
			elmPImgRight.style.backgroundImage = "url(./shared/Images/TrialCounter/pb_endred.gif)";
		}	
		else if ( nDaysRemain <= 14 )
		{
			strColourRGB = strOrangeRGB;
			// orange progressbar
			elmTopSlider.style.backgroundImage = "url(./shared/Images/TrialCounter/pborange.gif)";
			elmPImgLeft.style.backgroundImage = "url(./shared/Images/TrialCounter/pb_startorange.gif)";
			elmPImgRight.style.backgroundImage = "url(./shared/Images/TrialCounter/pb_endorange.gif)";
			
		}
		else
		{
			strColourRGB = strGrnRGB;
			// green progressbar
			elmTopSlider.style.backgroundImage = "url(./shared/Images/TrialCounter/pb.gif)";
			elmPImgLeft.style.backgroundImage = "url(./shared/Images/TrialCounter/pb_start.gif)";
			elmPImgRight.style.backgroundImage = "url(./shared/Images/TrialCounter/pb_end.gif)";
		}
		// change the colour of the numerical text
		if(null != elmCounter)
		{
			elmCounter.style.color = strColourRGB;
		}

		if(null != elmTopSlider)
		{
			// get the days remaining as a percentage
			var nPercent = (nDaysRemain / nTrialTotal);
			if(nPercent > 1)
			{
				nPercent = 1;
			}
			
			// get the total length of the progress bar and computer the desired
			// length of the slider
			if(null != elmBotSlider)
			{
				var nTotalWidth = parseInt(elmBotSlider.offsetWidth);
				var nPartWidth = (nTotalWidth * nPercent);
				if ( nPartWidth < nTotalWidth )
				{
					elmTopSlider.style.width = String(nPartWidth);
				}
			}
		}
	}
	else
	{
	    
		var tdPromo = document.getElementById("tdPromo")
		if (tdPromo != null)
		{
			tdPromo.style.display = "none";
		}
		// hide counter
		var divTrialRemain = document.getElementById("divTrialRemainCounter")
		if (divTrialRemain != null)
		{
			divTrialRemain.style.display = "none";
		}
		// hide progressbar
		var tdProgress = document.getElementById("tdProgress")
		if (tdProgress != null)
		{
			tdProgress.style.display = "none";
		}	
	}
	
	// Hide regSection
	var strURLlow = strURL.toLowerCase();
	var regType = extractParam(strURLlow, "regtype");
	var registered = extractParam(strURLlow, "isregistered");		
	var nStatus = extractParam(strURLlow, "status");
	var regSection = document.getElementById("regSection");
	if ( regSection != null )
	{
		var bHide = true;


			
		if ( regType != null && regType == "2" && registered != null && registered != "1" )
		{
			if (nStatus != null && nStatus != '3')
			{
				bHide = false;
			}	
		}
		
		if (bHide)
		{
			regSection.style.visibility = "hidden";
		}
		else
		{
			regSection.style.visibility = "visible";
		}	
	}
}

function ShowBanner(strShow)
{
	var banner = document.getElementById("banner");
	if (banner != null) 
	{
		if ( strShow.toLowerCase() == "no" )
		{
			banner.style.display = "none";
		}
		else 
		{
			banner.style.display = "inline";
		}
	}
}

function SetActivationCodeElement()
{
	var strActCode = getParamVal("activationcode");
	var oTxt = document.getElementById("txtActivationCode");
	
	if (strActCode != null && oTxt != null)
		oTxt.value = strActCode.toUpperCase();
}

function getParamString()
{
  var sParamString = "";
	var sQueryVar = location.search;
	var nIndex = sQueryVar.indexOf("?");
	if ( nIndex >= 0 )
	{
	  sParamString = sQueryVar.substr(nIndex);
  }	
  
  return sParamString;
}

function getParamVal(paramName)
{
	var sQueryVar = location.search.toLowerCase();
	sQueryVar = unescape(sQueryVar);
	return extractParam(sQueryVar, paramName);
}	

function extractParam(sQueryVar, paramName)
{
	if (sQueryVar == null || paramName == null)
	  return "";

	var nIndex = sQueryVar.indexOf("?");
	if (nIndex >= 0)
		sQueryVar = sQueryVar.substr(nIndex + 1);

	paramName = paramName.toLowerCase();

	var asParams = sQueryVar.split("&");
	var i;

	for (i in asParams)
	{
		var pair = asParams[i].split("=");
		if (pair.length == 2)
		{
			if (pair[0].toLowerCase() == paramName)
			{
				return pair[1];
			}
		}
	}

	return "";
}



function getParamValCaseSentitive(paramName)
{
	var sQueryVar = location.search;
	sQueryVar = unescape(sQueryVar);
	return extractParam(sQueryVar, paramName);
}


function EnableContinue(bEnable)
{
	var oBtn = document.getElementById("pcubtncontinue");
	if (oBtn != null)
	{
		oBtn.disabled = (bEnable == false);
	}
}

function ExpandRelativeURL(sPage)
{
	var sPath = location.pathname;
	var lastbkslsh = sPath.lastIndexOf("\\");
	if (lastbkslsh < 0)
	{
		lastbkslsh = sPath.lastIndexOf("/");
	}
	if (lastbkslsh < 0)
	{
		lastbkslsh = 0;
	}
	sPath = sPath.substr(0, lastbkslsh + 1);
	var sReturn = location.protocol + "//" + location.host + sPath + sPage;
	return sReturn;
}

function ValidateSerial(sSerial)
{
	if (sSerial != null && sSerial.match(/^\s*[a-zA-Z0-9]{2}\d{2}[a-zA-Z]{1}[a-zA-Z0-9]{2}-?[a-zA-Z0-9]{7}-?[a-zA-Z0-9]{7}-?[a-zA-Z0-9]{7}(?:\-?[a-zA-Z0-9]{5})?\s*$/) != null)
	{
		return true;
	}
	return false;
}

function ValidateActCode(sCode)
{
	if (sCode != null && sCode.match(/^\s*(?:[a-fA-F0-9]{4}-?){5}\s*$/) != null)
	{
		return true;
	}
	
	return false;
}

function FormatSerial(sCode)
{
	if (sCode != null && sCode.length >= 17 && sCode.indexOf("-") == -1)
	{
		var sNewCode = sCode.substr(0,7) + "-" + sCode.substr(7,7) + "-" + sCode.substr(14);
		sNewCode = sNewCode.toUpperCase();
		return sNewCode;
	}
	else
	{
		return sCode;
	}
}

function FormatActivationCode(sCode)
{
	if (sCode != null && sCode.length == 20 && sCode.indexOf("-") == -1)
	{
		var sNewCode = sCode.substr(0,4) + "-" + sCode.substr(4,4) + "-" + 
									 sCode.substr(8,4) + "-" + sCode.substr(12,4) + "-" + sCode.substr(16,4);
		sNewCode = sNewCode.toUpperCase();
		return sNewCode;
	}
	else
	{
		return sCode;
	}
}


function regChangeLaunchButtonTitle()
{
	var btnClose = document.getElementById("btnClose");
	if ( btnClose != null )
	{
		btnClose.style.display = "inline";
		
		var btnLaunch = document.getElementById("btnLaunch");
		if ( btnLaunch != null )
		{
			btnLaunch.style.display = "none";
		}
	}
}

function OnUnlock()
{
	location='pcucmd://Unlock?admin=yes'; 
	location='pcucmd://Next?condition=true'			
}

// Flash Player Version Detection - Rev 1.6
// Detect Client Browser type
// Copyright(c) 2005-2006 Adobe Macromedia Software, LLC. All rights reserved.
var isIE  = (navigator.appVersion.indexOf("MSIE") != -1) ? true : false;
var isWin = (navigator.appVersion.toLowerCase().indexOf("win") != -1) ? true : false;
var isOpera = (navigator.userAgent.indexOf("Opera") != -1) ? true : false;

function ControlVersion()
{
	var version;
	var axo;
	var e;

	// NOTE : new ActiveXObject(strFoo) throws an exception if strFoo isn't in the registry

	try {
		// version will be set for 7.X or greater players
		axo = new ActiveXObject("ShockwaveFlash.ShockwaveFlash.7");
		version = axo.GetVariable("$version");
	} catch (e) {
	}

	if (!version)
	{
		try {
			// version will be set for 6.X players only
			axo = new ActiveXObject("ShockwaveFlash.ShockwaveFlash.6");

			// installed player is some revision of 6.0
			// GetVariable("$version") crashes for versions 6.0.22 through 6.0.29,
			// so we have to be careful.

			// default to the first public version
			version = "WIN 6,0,21,0";

			// throws if AllowScripAccess does not exist (introduced in 6.0r47)
			axo.AllowScriptAccess = "always";

			// safe to call for 6.0r47 or greater
			version = axo.GetVariable("$version");

		} catch (e) {
		}
	}

	if (!version)
	{
		try {
			// version will be set for 4.X or 5.X player
			axo = new ActiveXObject("ShockwaveFlash.ShockwaveFlash.3");
			version = axo.GetVariable("$version");
		} catch (e) {
		}
	}

	if (!version)
	{
		try {
			// version will be set for 3.X player
			axo = new ActiveXObject("ShockwaveFlash.ShockwaveFlash.3");
			version = "WIN 3,0,18,0";
		} catch (e) {
		}
	}

	if (!version)
	{
		try {
			// version will be set for 2.X player
			axo = new ActiveXObject("ShockwaveFlash.ShockwaveFlash");
			version = "WIN 2,0,0,11";
		} catch (e) {
			version = -1;
		}
	}

	return version;
}

// JavaScript helper required to detect Flash Player PlugIn version information
function GetSwfVer(){
	// NS/Opera version >= 3 check for Flash plugin in plugin array
	var flashVer = -1;

	if (navigator.plugins != null && navigator.plugins.length > 0) {
		if (navigator.plugins["Shockwave Flash 2.0"] || navigator.plugins["Shockwave Flash"]) {
			var swVer2 = navigator.plugins["Shockwave Flash 2.0"] ? " 2.0" : "";
			var flashDescription = navigator.plugins["Shockwave Flash" + swVer2].description;
			var descArray = flashDescription.split(" ");
			var tempArrayMajor = descArray[2].split(".");
			var versionMajor = tempArrayMajor[0];
			var versionMinor = tempArrayMajor[1];
			var versionRevision = descArray[3];
			if (versionRevision == "") {
				versionRevision = descArray[4];
			}
			if (versionRevision[0] == "d") {
				versionRevision = versionRevision.substring(1);
			} else if (versionRevision[0] == "r") {
				versionRevision = versionRevision.substring(1);
				if (versionRevision.indexOf("d") > 0) {
					versionRevision = versionRevision.substring(0, versionRevision.indexOf("d"));
				}
			}
			var flashVer = versionMajor + "." + versionMinor + "." + versionRevision;
		}
	}
	// MSN/WebTV 2.6 supports Flash 4
	else if (navigator.userAgent.toLowerCase().indexOf("webtv/2.6") != -1) flashVer = 4;
	// WebTV 2.5 supports Flash 3
	else if (navigator.userAgent.toLowerCase().indexOf("webtv/2.5") != -1) flashVer = 3;
	// older WebTV supports Flash 2
	else if (navigator.userAgent.toLowerCase().indexOf("webtv") != -1) flashVer = 2;
	else if ( isIE && isWin && !isOpera ) {
		flashVer = ControlVersion();
	}
	return flashVer;
}

// When called with reqMajorVer, reqMinorVer, reqRevision returns true if that version or greater is available
function DetectFlashVer(reqMajorVer, reqMinorVer, reqRevision)
{
	versionStr = GetSwfVer();
	if (versionStr == -1 ) {
		return false;
	} else if (versionStr != 0) {
		if(isIE && isWin && !isOpera) {
			// Given "WIN 2,0,0,11"
			tempArray         = versionStr.split(" "); 	// ["WIN", "2,0,0,11"]
			tempString        = tempArray[1];			// "2,0,0,11"
			versionArray      = tempString.split(",");	// ['2', '0', '0', '11']
		} else {
			versionArray      = versionStr.split(".");
		}
		var versionMajor      = versionArray[0];
		var versionMinor      = versionArray[1];
		var versionRevision   = versionArray[2];

        	// is the major.revision >= requested major.revision AND the minor version >= requested minor
		if (versionMajor > parseFloat(reqMajorVer)) {
			return true;
		} else if (versionMajor == parseFloat(reqMajorVer)) {
			if (versionMinor > parseFloat(reqMinorVer))
				return true;
			else if (versionMinor == parseFloat(reqMinorVer)) {
				if (versionRevision >= parseFloat(reqRevision))
					return true;
			}
		}
		return false;
	}
}

function AC_AddExtension(src, ext)
{
  if (src.indexOf('?') != -1)
    return src.replace(/\?/, ext+'?');
  else
    return src + ext;
}

function AC_Generateobj(objAttrs, params, embedAttrs)
{
    var str = '';
    if (isIE && isWin && !isOpera)
    {
  		str += '<object ';
  		for (var i in objAttrs)
  			str += i + '="' + objAttrs[i] + '" ';
  		for (var i in params)
  			str += '><param name="' + i + '" value="' + params[i] + '" /> ';
  		str += '></object>';
    } else {
  		str += '<embed ';
  		for (var i in embedAttrs)
  			str += i + '="' + embedAttrs[i] + '" ';
  		str += '> </embed>';
    }

    document.write(str);
}

function AC_FL_RunContent(){
  var ret =
    AC_GetArgs
    (  arguments, ".swf", "movie", "clsid:d27cdb6e-ae6d-11cf-96b8-444553540000"
     , "application/x-shockwave-flash"
    );
  AC_Generateobj(ret.objAttrs, ret.params, ret.embedAttrs);
}

function AC_GetArgs(args, ext, srcParamName, classid, mimeType){
  var ret = new Object();
  ret.embedAttrs = new Object();
  ret.params = new Object();
  ret.objAttrs = new Object();
  for (var i=0; i < args.length; i=i+2){
    var currArg = args[i].toLowerCase();

    switch (currArg){
      case "classid":
        break;
      case "pluginspage":
        ret.embedAttrs[args[i]] = args[i+1];
        break;
      case "src":
      case "movie":
        args[i+1] = AC_AddExtension(args[i+1], ext);
        ret.embedAttrs["src"] = args[i+1];
        ret.params[srcParamName] = args[i+1];
        break;
      case "onafterupdate":
      case "onbeforeupdate":
      case "onblur":
      case "oncellchange":
      case "onclick":
      case "ondblClick":
      case "ondrag":
      case "ondragend":
      case "ondragenter":
      case "ondragleave":
      case "ondragover":
      case "ondrop":
      case "onfinish":
      case "onfocus":
      case "onhelp":
      case "onmousedown":
      case "onmouseup":
      case "onmouseover":
      case "onmousemove":
      case "onmouseout":
      case "onkeypress":
      case "onkeydown":
      case "onkeyup":
      case "onload":
      case "onlosecapture":
      case "onpropertychange":
      case "onreadystatechange":
      case "onrowsdelete":
      case "onrowenter":
      case "onrowexit":
      case "onrowsinserted":
      case "onstart":
      case "onscroll":
      case "onbeforeeditfocus":
      case "onactivate":
      case "onbeforedeactivate":
      case "ondeactivate":
      case "type":
      case "codebase":
        ret.objAttrs[args[i]] = args[i+1];
        break;
      case "id":
      case "width":
      case "height":
      case "align":
      case "vspace":
      case "hspace":
      case "class":
      case "title":
      case "accesskey":
      case "name":
      case "tabindex":
        ret.embedAttrs[args[i]] = ret.objAttrs[args[i]] = args[i+1];
        break;
      default:
        ret.embedAttrs[args[i]] = ret.params[args[i]] = args[i+1];
    }
  }
  ret.objAttrs["classid"] = classid;
  if (mimeType) ret.embedAttrs["type"] = mimeType;
  return ret;
}

/* function used for tab switching */

function switchtabs(id,exc,max)
{
	//exception where 'exc' tab doesnt show bottom table
	if (exc > 0)
	{
		var bottomtable = document.getElementById('bottomtable');
		
		if (id == exc) {
			bottomtable.style.display = 'none';
		}
		else {
			bottomtable.style.display = '';
		}
	}

	for(var i=1;i<=max;i++)
	{
		var tabcontent = document.getElementById('tabpage' + i);
		var tabitem = document.getElementById('tab' + i);

		//change tab content
		if (i!=id) {
			tabcontent.style.display = 'none';
		}
		else {
			tabcontent.style.display = '';
		}

		//highlight tab
		if (i!=id) {
			tabitem.style.backgroundColor = 'black';
			tabitem.childNodes[0].style.color = 'white';
		}
		else {
			tabitem.style.backgroundColor = 'white';
			tabitem.childNodes[0].style.color = 'black';
		}
	}
}

/* function used for opening up a URL or predefined content in a pop-up window */

function popitup(url,height,width)
{	
	var dimensions = "height=" + height + ",width=" + width + "'";

	//generate new window with a unique name everytime someone clicks
	var uniqid = Math.floor(Math.random()*100001);
	newwindow = window.open(url,'popup' + uniqid,dimensions);

	//popups with predefined content
	if (url == 'about:blank')
	{
		if (arguments[3] == 'MB') {
			newwindow.document.body.innerHTML = '<h2 style="font-size: 14px">Corel 30-Day Money-Back Guarantee</h2><div>At Corel, we\'re proud of our products and we\'re confident that you\'ll be satisfied with your purchase. That\'s why we offer the <strong>Corel 30-Day Money-Back Guarantee</strong>. Buy today and if you\'re not completely satisfied with your purchase, Corel will refund your money based on your original method of payment &#151; excluding any shipping costs. When you contact us about a return, your Corel Customer Service Representative may inquire as to the reason for your return so we can improve our products and our customer experience.</div><div style="margin-top: 10px"><strong>Please note</strong>: Only products purchased directly from Corel qualify for the Corel 30-Day Money-Back Guarantee. Purchases made from a retailer must be returned to where you made your purchase and are subject to the return policy of that retailer. Purchases of transactional license, maintenance purchases, Wilcom DecoStudio products, Wacom Cintiq products and download insurance are excluded from the 30 day money back guarantee.</div>';
		}
		else if (arguments[3] == 'UPGRADE') {
			newwindow.document.body.innerHTML = '<h2 style="font-size: 14px">Here\'s how you qualify to save big!</h2><div>If you own any of the following products, you can get Corel Paint Shop Pro Photo X2 Ultimate for only $59.99:</div><ul><li>Microsoft&reg; Digital Image Pro, Digital Image Suite or Picture It</li><li>Ulead&reg;/Nova&reg; Photo Impact or Nova&reg; Photo Explosion</li><li>CorelDRAW&reg;, CorelDRAW&reg; Graphics Suite, Corel&reg; Painter or Corel&reg; Painter Essentials</li></ul><div style="font-weight: bold">Get upgrade savings and order Paint Shop Pro Photo X2 Ultimate today!</div>';
		}
		else {
			newwindow.document.body.innerHTML = arguments[3];
		}
	}
	newwindow.focus();
}

function setCookie(name, value) 
{
    var newCookie = name + "=" + escape(value) + "; path=/";
    document.cookie = newCookie;
}

function getCookie(name) 
{
    var prefix = name + "="
    var cookieStartIndex = document.cookie.indexOf(prefix)
    if (cookieStartIndex == -1)
        return null
    var cookieEndIndex = document.cookie.indexOf(";", cookieStartIndex + prefix.length)
    if (cookieEndIndex == -1)
        cookieEndIndex = document.cookie.length
    return unescape(document.cookie.substring(cookieStartIndex + prefix.length, cookieEndIndex))
}

function deleteCookie(name) 
{
    if (getCookie(name)) 
    {
        document.cookie = name + "=; path=/; expires=Thu, 01-Jan-70 00:00:01 GMT";
    }
}

function removeAllSpaces(strVal)
{
    return strVal.replace(/\s/g,"");
}

function init_trial_page()
{
	try {
		SetExtensionPeriodText();
		headInit();
	}
	catch(e) {};
}

function SetExtensionPeriodText()
{
	var extensionPeriodElmentList = document.getElementsByName("ExtensionPeriodText");
	if ( extensionPeriodElmentList != null && extensionPeriodElmentList.length > 0)
	{
		var extensionPeriodDays = getParamVal("ExtensionPeriod");
		if ( extensionPeriodDays!= null)
		{
			var extensionPeriodElment = null;
			var i = 0;
			for ( i = 0; i < extensionPeriodElmentList.length; i++)
			{
				extensionPeriodElment = extensionPeriodElmentList[i];
				if ( extensionPeriodElment != null)
				{
					extensionPeriodElment.innerText = extensionPeriodDays;
				}
			}
		}
	}
}

/* the caller to the below 2 functions has to implement the callback functions calBack(strReturn) which will be called upon completion of the async call */
function sendEmailValidation() {
    var xmlhttp;
	//http://www.corel.com/corel/ipm/sendEmailValidation.jsp was ATG url
    var strURL = "http://iws.corel.com/ipmws/sendEmailValidation" + getParamString();
    if (window.XMLHttpRequest) {
        // for IE7 and higher
        xmlhttp = new XMLHttpRequest();
    }
    else {
        // for IE6
        xmlhttp = new ActiveXObject("Microsoft.XMLHTTP");
    }
    xmlhttp.onreadystatechange = function () {
        var strReturn = "";
        if (xmlhttp.readyState == 4) {
            if (xmlhttp.status == 200) {
                strReturn = xmlhttp.responseText;
                var token = "<body>";
                var tokenIndex = strReturn.indexOf(token);
                if (tokenIndex > -1) {
                    strReturn = strReturn.substring(tokenIndex + token.length);
                    token = "</body>";
                    tokenIndex = strReturn.indexOf(token);
                    if (tokenIndex > -1) {
                        strReturn = removeAllSpaces(strReturn.substring(0, tokenIndex));
                        callBack(strReturn);
                    }
                }
            }
            else {
                callBack("1");    //failed to send/connect
            }
        }
    }
    xmlhttp.open("GET", strURL, true);
    xmlhttp.send();
}

function sendPremTrialStart() {
    var xmlhttp;
	var strURL = "https://iws.corel.com/ipmws/redirect/index.jsp" + getParamString();
    if (window.XMLHttpRequest) {
        // for IE7 and higher
        xmlhttp = new XMLHttpRequest();
    }
    else {
        // for IE6
        xmlhttp = new ActiveXObject("Microsoft.XMLHTTP");
    }
    xmlhttp.onreadystatechange = function () {
        var strReturn = "";
        if (xmlhttp.readyState == 4) {
            if (xmlhttp.status == 200) {
                strReturn = xmlhttp.responseText;
                var token = "<body>";
                var tokenIndex = strReturn.indexOf(token);
                if (tokenIndex > -1) {
                    strReturn = strReturn.substring(tokenIndex + token.length);
                    token = "</body>";
                    tokenIndex = strReturn.indexOf(token);
                    if (tokenIndex > -1) {
                        strReturn = removeAllSpaces(strReturn.substring(0, tokenIndex));
                        callBack(strReturn);
                    }
                }
            }
            else {
                callBack("1");    //failed to send/connect
            }
        }
    }
    xmlhttp.open("GET", strURL, true);
    xmlhttp.send();
}

/* email validation reminders */

function changeEmail() {
	//https://www.corel.com/corel/accountEditEmail.jsp was ATG URL
    location = "pcucmd://fullbrowse?url='https://account.corel.com/accounts/user/accountEditEmail'";
}

/* form label functions */
function fieldin(fieldobj) {
	fieldobj.style.backgroundColor = '#fcfcfc';
}

function fieldout(fieldobj)
{
	if(fieldobj.value == '') {
		fieldobj.style.backgroundColor = 'transparent';
	}
}
/* end form label functions */

function getCorelComServer()
{
	var s = getParamVal("CorelComServer");
	if (s == '') {
		s = "http://www.corel.com/";
	}
	return s;
}

function launchStoreEx(noAdmin, fullBrowser, w, h) {
	var storeUrl = addGoogleAnalyticsParamsToUrl("http://apps.corel.com/redirect?_redirect=pcupurchase" + getPCUNVPs(true, noAdmin));
	var url = "pcucmd://browse?URL='" + storeUrl + "'&width=" + w + "&height=" + h;
	if(fullBrowser == "true" || fullBrowser == true){
			url = "pcucmd://fullbrowse?url='" + storeUrl + "'";
	}
	document.location.href = url;
}

// start moving window
function StartMoveWindow() {
   location='pcucmd://StartRCMoveWindow/';
}

//check if page inside iframe is on a different domain
function canAccessIframe(iframe) {
	try {
		var doc = iframe.contentWindow.document;
    }
	catch(err){
      return false;
    }
	return true;
}

//close container dialog
function CloseDlg() {

    //first check to see if pcubtncancel is in the button bar if the container
    var cancelbtn = document.getElementById("pcubtncancel");
    if (cancelbtn != null && cancelbtn.onclick) {
        cancelbtn.onclick({ target: cancelbtn });
    }
    else {
        //check if there is a content iframe
        var ifrm = document.getElementById("Content");
        if (ifrm != null) {
            //check if the content iframe document can be accessed
            if (canAccessIframe(ifrm) == true) {
                //check if content iframe has a cancel button with a custom message
                cancelbtn = ifrm.contentWindow.document.getElementById("pcubtncancel");
                if (cancelbtn != null && cancelbtn.onclick) {
                    cancelbtn.onclick({ target: cancelbtn });
                }
                else {
                    document.location = 'pcucmd://cancel';
                }
            }
            else {
                document.location = 'pcucmd://cancel';
            }
        }
        else {
            document.location = 'pcucmd://cancel';
        }

    } //document.getElementById("pcubtncancel"); - button bar
}

//onclick forgot password
function forgotPassword() {
	location.href = "pcucmd://Connect";
	location.href = "pcucmd://setsequence?name=forgot_password&condition=true";
}

function supportlink() {
	location.href = "pcucmd://fullbrowse?url='<CorelComServer>/membersupport'";
	return false;
}
