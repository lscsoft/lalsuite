<!-- start_time -->
<div id="time">
<input type="radio" name="starttimevsgps" value="gps">
Start GPS: <input type="text" size=10 maxlength=10 name="startgps"> or 
<input type="radio" name="starttimevsgps" value="time" checked="checked">
            Start time: <select name="startmonth"><?php require '/usr1/ldbd/glue/var/php/script/form_month_list.php' ?></select>
                     <select name="startday"><?php require '/usr1/ldbd/glue/var/php/script/form_day_list.php' ?></select>
                     <select name="startyear"><?php require '/usr1/ldbd/glue/var/php/script/form_year_list.php' ?></select>
            <input type="text" name="starthour" value=00 size=2 maxlength=2>
            :<input type="text" name="startmin" value=00 size=2 maxlength=2>
            :<input type="text" name="startsec" value=00 size=2 maxlength=2>
            <select name="starttype">
                    <option value="24Hour">24 Hr</option>
                    <option value="am">am</option>
                    <option value="pm">pm</option>
            </select>
            <select name="startzone">
                    <option value="UTC">UTC</option>
                    <option value="Central">Central</option>
                    <option value="Pacific">Pacific</option>
            </select>
<div class="placeholder_bg"></div>

<!-- end_time -->
<input type="radio" name="stoptimevsgps" value="gps">
Stop GPS: <input type="text" size=10 maxlength=10 name="stopgps"> or
<input type="radio" name="stoptimevsgps" value="time" checked="checked">
           Stop time: <select name="stopmonth"><?php require '/usr1/ldbd/glue/var/php/script/form_month_list.php' ?></select>
                <select name="stopday"><?php require '/usr1/ldbd/glue/var/php/script/form_day_list.php' ?></select>
                <select name="stopyear"><?php require '/usr1/ldbd/glue/var/php/script/form_year_list.php' ?></select>
        <input type="text" name="stophour" value=00 size=2 maxlength=2>
        :<input type="text" name="stopmin"  value=00 size=2 maxlength=2>
        :<input type="text" name="stopsec"  value=00 size=2 maxlength=2>
        <select name="stoptype">
                    <option value="24Hour">24 Hr</option>
                    <option value="am">am</option>
                    <option value="pm">pm</option>
        </select>
        <select name="stopzone">
                    <option value="UTC">UTC</option>
                    <option value="Central">Central</option>
                    <option value="Pacific">Pacific</option>
        </select>

</div>
